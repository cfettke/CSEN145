#include <vector>
#include <fstream>

class Vector
{
private:
    
public:
    int n;
    std::vector<double> buffer;
    Vector(int n_values=100);
    ~Vector();
};

// Initialize vector with random values
Vector::Vector(int n_values)
{

    srand(1337);
    n = n_values;
    for(int i=0; i<n; ++i){
        float value = static_cast <float> (rand()) / static_cast <float> (rand());
        this->buffer.push_back(value);
    }

    // Write to File
    std::ofstream fp("vec.txt");
    for(auto &i : buffer){
        fp << i << std::endl;
    }
    fp.close();

}

Vector::~Vector()
{
}
