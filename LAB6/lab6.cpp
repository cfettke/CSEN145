#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cstring>      /* strcasecmp */
#include <cstdint>
#include <assert.h>
#include <vector>       // std::vector
#include <algorithm>    // std::random_shuffle
#include <random>
#include <stdexcept>
#include <mpi.h>
#include <cmath>
#include <sstream>
#include <limits>

using namespace std;

using idx_t = std::uint32_t;
using val_t = float;
using ptr_t = std::uintptr_t;

/**
 * CSR structure to store search results
 */
typedef struct csr_t {
  idx_t nrows; // number of rows
  idx_t ncols; // number of rows
  idx_t * ind; // column ids
  val_t * val; // values
  ptr_t * ptr; // pointers (start of row in ind/val)

  csr_t()
  {
    nrows = ncols = 0;
    ind = nullptr;
    val = nullptr;
    ptr = nullptr;
  }

  /**
   * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
   * @param nrows Number of rows
   * @param nnz   Number of non-zeros
   */
  void reserve(const idx_t nrows, const ptr_t nnz)
  {
    if(nrows > this->nrows){
      if(ptr){
        ptr = (ptr_t*) realloc(ptr, sizeof(ptr_t) * (nrows+1));
      } else {
        ptr = (ptr_t*) malloc(sizeof(ptr_t) * (nrows+1));
        ptr[0] = 0;
      }
      if(!ptr){
        throw std::runtime_error("Could not allocate ptr array.");
      }
    }
    if(ind){
      ind = (idx_t*) realloc(ind, sizeof(idx_t) * nnz);
    } else {
      ind = (idx_t*) malloc(sizeof(idx_t) * nnz);
    }
    if(!ind){
      throw std::runtime_error("Could not allocate ind array.");
    }
    if(val){
      val = (val_t*) realloc(val, sizeof(val_t) * nnz);
    } else {
      val = (val_t*) malloc(sizeof(val_t) * nnz);
    }
    if(!val){
      throw std::runtime_error("Could not allocate val array.");
    }
    this->nrows = nrows;
  }

  csr_t ( const csr_t &other)
  {
    this->nrows = this->ncols = 0;
    this->ptr = nullptr;
    this->ind = nullptr;
    this->val = nullptr;
    this->reserve(other.nrows, other.ptr[other.nrows]);
    memcpy(ptr, other.ptr, sizeof(ptr_t) * (nrows+1));
    memcpy(ind, other.ind, sizeof(idx_t) * ptr[nrows]);
    memcpy(val, other.val, sizeof(val_t) * ptr[nrows]);
    this->ncols = other.ncols;
  }

  /**
   * Create random matrix with given sparsity factor.
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param factor   Sparsity factor
   */
  static csr_t * random(const idx_t nrows, const idx_t ncols, const double factor)
  {
    ptr_t nnz = (ptr_t) (factor * nrows * ncols);
    if(nnz >= nrows * ncols / 2.0){
      throw std::runtime_error("Asking for too many non-zeros. Matrix is not sparse.");
    }
    auto mat = new csr_t();
    mat->reserve(nrows, nnz);
    mat->ncols = ncols;

    /* fill in ptr array; generate random row sizes */
    unsigned int seed = (unsigned long) mat;
    long double sum = 0;
    for(idx_t i=1; i <= mat->nrows; ++i){
      mat->ptr[i] = rand_r(&seed) % ncols;
      sum += mat->ptr[i];
    }
    for(idx_t i=0; i < mat->nrows; ++i){
      double percent = mat->ptr[i+1] / sum;
      mat->ptr[i+1] = mat->ptr[i] + (ptr_t)(percent * nnz);
      if(mat->ptr[i+1] > nnz){
        mat->ptr[i+1] = nnz;
      }
    }
    if(nnz - mat->ptr[mat->nrows-1] <= ncols){
      mat->ptr[mat->nrows] = nnz;
    }

    /* fill in indices and values with random numbers */
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      unsigned int seed = (unsigned long) mat * (1+tid);
      std::vector<int> perm;
      for(idx_t i=0; i < ncols; ++i){
        perm.push_back(i);
      }
      std::random_device seeder;
      std::mt19937 engine(seeder());

      #pragma omp for
      for(idx_t i=0; i < nrows; ++i){
        std::shuffle(perm.begin(), perm.end(), engine);
        for(ptr_t j=mat->ptr[i]; j < mat->ptr[i+1]; ++j){
          mat->ind[j] = perm[j - mat->ptr[i]];
          mat->val[j] = ((double) rand_r(&seed)/rand_r(&seed));
        }
      }
    }

    return mat;
  }

  string info(const string name="") const {
    return (name.empty() ? "CSR" : name) + "<" + to_string(nrows) + ", " + to_string(ncols) + ", " +(ptr ? to_string(ptr[nrows]) : "0") + ">";
  }

  /** 
   * Read the matrix from a CLUTO file.
   * The first line is "nrows ncols nnz".
   * Each other line is a row in the matrix containing column ID-value pairs for non-zeros in the row.
  */
  void read(const std::string &filename)
  {
    FILE * infile = fopen(filename.c_str(), "r");
    char * line = NULL;
    size_t n, nr, nnz;
    char *head;
    char *tail;
    idx_t cid;
    double dval;
    
    if (!infile) {
      throw std::runtime_error("Could not open CLU file\n");
    }
    if(getline (&line, &n, infile) < 0){
      throw std::runtime_error("Could not read first line from CLU file\n");
    }
    //read matriz size info
    size_t rnrows, rncols, rnnz;
    sscanf(line, "%zu %zu %zu", &rnrows, &rncols, &rnnz);

    //allocate space
    this->reserve(rnrows, rnnz);
    ncols = rncols;
    
    //read in rowval, rowind, rowptr
    this->ptr[0]= 0;
    nnz = 0;
    nr = 0;

    while(getline(&line, &n, infile) != -1){
      head = line;
      while (1) {
        cid = (idx_t) strtol(head, &tail, 0);
        if (tail == head)
          break;
        head = tail;

        if(cid <= 0){
          throw std::runtime_error("Invalid column ID while reading CLUTO matrix\n");
        }
        this->ind[nnz] = cid - 1; //csr/clu files have 1-index based column IDs and our matrix is 0-based.
        dval = strtod(head, &tail);
        head = tail;
        this->val[nnz++] = dval;
      }
      this->ptr[nr+1] = nnz;
      nr++;
    }
    assert(nr == rnrows);
    free(line);
    fclose(infile);
  }
  
  /** 
   * Read the matrix from a CLUTO file.
   * The first line is "nrows ncols nnz".
   * Each other line is a row in the matrix containing column ID-value pairs for non-zeros in the row.
  */
  static csr_t * from_CLUTO(const std::string &filename)
  {
    auto mat = new csr_t();
    mat->read(filename);
    return mat;
  }

  /**
   * Write matrix to text file
   * @param output_fpath File to write to
   */
  void write(const std::string output_fpath, const bool header=false)
  {
    std::fstream resfile;
    resfile.open(output_fpath, std::ios::out);
    if(!resfile){
      throw std::runtime_error("Could not open output file for writing.");
    }
    if(header){
      resfile << nrows << " " << ncols << " " << ptr[nrows] << std::endl;
    }
    for(idx_t i=0; i < nrows; ++i){
      for(ptr_t j=ptr[i]; j < ptr[i+1]; ++j){
        resfile << ind[j] << " " << val[j];
        if(j+1 < ptr[i+1]){
          resfile << " ";
        }
      }
      //resfile << std::endl;
    }
    resfile.close();
  }

  ~csr_t()
  {
    if(ind){
      free(ind);
    }
    if(val){
      free(val);
    }
    if(ptr){
      free(ptr);
    }
  }
} csr_t;

/**
 * Ensure the matrix is valid
 * @param mat Matrix to test
 */
void test_matrix(csr_t * mat){
  auto nrows = mat->nrows;
  auto ncols = mat->ncols;
  assert(mat->ptr);
  auto nnz = mat->ptr[nrows];
  for(idx_t i=0; i < nrows; ++i){
    assert(mat->ptr[i] <= nnz);
  }
  for(ptr_t j=0; j < nnz; ++j){
    assert(mat->ind[j] < ncols);
  }
}

void writeClutoFile(const std::vector<idx_t>& row_ptr, const std::vector<idx_t>& col_idx, const std::vector<val_t>& values, const std::string& filename) {
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Write the header section
    int nrows = row_ptr.size() - 1;
    int nnz = row_ptr[nrows];
    int ncols = 0;  // Calculate the actual number of columns

    for (idx_t i = 0; i < nnz; ++i) {
        if (col_idx[i] > ncols) {
            ncols = col_idx[i];
        }
    }

    file << nrows << " " << ncols << " " << nnz << endl;

    // Write the data section
    for (idx_t i = 0; i < nrows; ++i) {
        for (idx_t j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            file << col_idx[j] << " " << values[j] << " ";
        }
        file << endl;
    }

    file.close();
}

/**
 * Sends ptr vector, index vector, and vals vector
 * @param rowPtr RowPtr vector to send
 * @param columnInd vector to send
 * @param rowVal vector to send
 * @param destinaiton desination process rank
*/
void sendVectors(const vector<idx_t>& rowPtr, const vector<idx_t>& columnInd, const vector<val_t>& rowVal,int destination) {
    MPI_Send((void*)rowPtr.data(), rowPtr.size() * sizeof(idx_t), MPI_BYTE, destination, 0, MPI_COMM_WORLD);
    MPI_Send((void*)columnInd.data(), columnInd.size() * sizeof(idx_t), MPI_BYTE, destination, 1, MPI_COMM_WORLD);
    MPI_Send((void*)rowVal.data(), rowVal.size() * sizeof(val_t), MPI_BYTE, destination, 2, MPI_COMM_WORLD);
}

/**
 * Recieves incoming vectors via mpi recieve
 * @param incomingRowPtr RowPtr vector to store received info
 * @param incomingColumnInd Column Indices Vector
 * @param incomingRowVal Values Vector
 * @param processRank Process Rank of caller process
*/
void receiveVectors(vector<idx_t>& incomingRowPtr, vector<idx_t>& incomingColumnInd, vector<val_t>& incomingRowVal, int processRank) {
    MPI_Status status;

    //gets values and MPI info
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    int incomingRowPtr_size;
    MPI_Get_count(&status, MPI_BYTE, &incomingRowPtr_size);
    incomingRowPtr.resize(incomingRowPtr_size / sizeof(idx_t));
    MPI_Recv((void*)incomingRowPtr.data(), incomingRowPtr_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    idx_t nrows = incomingRowPtr.back(); 
    incomingRowPtr.pop_back();  //used because number of ptrs is stored in last value

    MPI_Probe(0, 1, MPI_COMM_WORLD, &status);
    int incomingColumnInd_size;
    MPI_Get_count(&status, MPI_BYTE, &incomingColumnInd_size);
    incomingColumnInd.resize(incomingColumnInd_size / sizeof(idx_t));
    MPI_Recv((void*)incomingColumnInd.data(), incomingColumnInd_size, MPI_BYTE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    idx_t ncols = incomingColumnInd.back(); 
    incomingColumnInd.pop_back(); //used because number of inds is stored in last value

    //Receive incomingRowVal
    MPI_Probe(0, 2, MPI_COMM_WORLD, &status);
    int incomingRowVal_size;
    MPI_Get_count(&status, MPI_BYTE, &incomingRowVal_size);
    incomingRowVal.resize(incomingRowVal_size / sizeof(val_t));
    MPI_Recv((void*)incomingRowVal.data(), incomingRowVal_size, MPI_BYTE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //Copy data to csr
    string fileName = "process" + to_string(processRank) + "mpi.clu";
    writeClutoFile(incomingRowPtr, incomingColumnInd, incomingRowVal, fileName);
    auto processCSR = csr_t::from_CLUTO(fileName);
    test_matrix(processCSR);
    cout << processCSR->info() << endl;

    delete processCSR;
}

/**
 * Partitions cluto into evenly distributed number of elemnts based on number of processes
 * @param filename name of cluto file
 * @param numProcesses number of processes
 * @param processRank rank of each process
*/
void partition_CLUTO(char* filename, int numProcesses, int processRank) {    

    //root process divides matrix into chunks
    if(processRank == 0){
      ifstream file(filename);

      if (!file.is_open()) {
          cerr << "Error opening file: " << filename << endl;
          return;
      }

      //gets rows, cols, and non zeros from first line then reads first line
      idx_t nrows, ncols, nnz;
      file >> nrows >> ncols >> nnz;
      string line;
      getline(file, line);

      //evenly distributes non zeros based on number of processes
      size_t nonzeros_per_process = ceil(static_cast<double>(nnz) / static_cast<double>(numProcesses-1));
      std::cout << "Partioned to " << numProcesses - 1 << " parts, " << nonzeros_per_process << " nonzeros per process\n";

      //stores beginning and end of each process
      cout << "Available processes: " << numProcesses << endl;
      vector<idx_t> values_per_process; 
      values_per_process.push_back(0);
      for(int i = 0; i < numProcesses - 1; i++) {
          values_per_process.push_back(nonzeros_per_process * (i + 1) - 1);
      }
      values_per_process.push_back(nnz - 1);

      idx_t nonzeros_allocated = 0;
      int destination = 1;

      //Create vectors for sending data
      vector<idx_t> rowPtr;
      vector<idx_t> columnInd;
      vector<val_t> rowVal;
      rowPtr.push_back(0);

      for (idx_t row = 0; row < nrows; ++row) {
        //For new row, read the line, and then iterate thru values
        getline(file, line); 
        stringstream ss(line); 
        
        //temp variables for values
        idx_t colIndexTemp;
        val_t valTemp;

        //iterate thru pairs
        while (ss >> colIndexTemp >> valTemp) {
          //checks if enough nonzeros have been iterated through to send to a process
          if(nonzeros_allocated >= nonzeros_per_process){
            //reset count
            nonzeros_allocated = 0;

            //updates nonzeros
            nonzeros_per_process = values_per_process[destination] - values_per_process[destination - 1];

            //adds col and rows to back so that info can be recovered without sending another mpi send
            columnInd.push_back(ncols);
            rowPtr.push_back(nrows);
            cout << destination << endl;

            //sends rowPtrs, columnInds, rowVals, and destination then increments desinaiton process
            sendVectors(rowPtr, columnInd, rowVal, destination);
            destination++;

            //set vectors to zero so that for the next process the initial rows are kept
            for(size_t zeros = 0; zeros < rowPtr.size(); zeros++){
              rowPtr[zeros] = 0;
            }
            columnInd.clear();
            rowVal.clear();
          }

          columnInd.push_back(colIndexTemp);
          rowVal.push_back(valTemp);
          nonzeros_allocated++;
        }
        rowPtr.push_back(nonzeros_allocated);
        //go to next row
      }
      //sends final set of vectors
      sendVectors(rowPtr, columnInd, rowVal, destination);
      file.close();
    } else {
      vector<idx_t> incomingRowPtr;
      vector<idx_t> incomingColumnInd;
      vector<val_t> incomingRowVal;
      receiveVectors(incomingRowPtr, incomingColumnInd, incomingRowVal, processRank);
    }
}

int main(int argc, char *argv[])
{
  //these are all intializiation steps
  MPI_Init(NULL, NULL); // initialize MPI environment

  int world_size; // number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank; // the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  //function to split file
  partition_CLUTO(argv[1], world_size, world_rank);

  MPI_Finalize();

  return 0;
} 