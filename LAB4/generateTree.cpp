#include <iostream>
#include <random>
#include <ctime>
#include <omp.h>
#include <queue>
#include <fstream>
#include <math.h>
#include <ctime>

using namespace std;

struct treeNode {
    public:
        double val;
        treeNode* left;
        treeNode* right;
        treeNode() : val(0.0), left(nullptr), right(nullptr) {}
        treeNode(double x) : val(x), left(nullptr), right(nullptr) {}
        treeNode(double x, treeNode* l, treeNode* r) : val(x), left(l), right(r) {}
};

void preOrderTraversalOMP(treeNode* root, int level, double sum) {
    //if(root != nullptr) {
        double temp;
        //omp_set_num_threads(4);
        //cout << root->val << ", level = " << level << end;
        #pragma omp parallel sections shared(sum)
        {
            #pragma omp section 
            {  
                temp = cos(pow(root->val, 5));
                temp = pow(temp, 10);
                cout << temp;// << ", num threads = " << omp_get_num_threads();
                for(int i = 0; i < 300; i++) {
                    cout << "-";
                }
                cout << endl;
            }
            #pragma omp section
            {
                if(root->left != nullptr) {
                    preOrderTraversalOMP(root->left, level + 1, sum);
                }
            }
            #pragma omp section
            {
                if(root->right != nullptr) {
                    preOrderTraversalOMP(root->right, level + 1, sum);
                }
            }
        }
        
    //}
}

double sum_recursive(treeNode* current, int depth) {
    if (current == NULL) {
        return 0.0;
    }
    else {
        double sum;
        if(depth <= 1) {
            double left, right, sum;
            #pragma omp task shared(left)
            {
                left = sum_recursive(current->left, depth + 1);
            }
            #pragma omp task shared(right)
            {
                right = sum_recursive(current->right, depth + 1);
            }
            #pragma omp taskwait
            return sin(current->val) + left + right;
        }
        else {
            sum = (sin(current->val) + sum_recursive(current->left, depth + 1) + sum_recursive(current->right, depth + 1));
            return sum;
        }
    }
}







void preOrderTraversal(treeNode* root, int level, double sum) {
    //if(root != nullptr) {
        double temp;
        //omp_set_num_threads(4);
        //cout << root->val << ", level = " << level << end;
        sum += root->val;
        temp = pow(root->val, 5);
        temp = pow(temp, 10);
        cout << temp;// << ", num threads = " << omp_get_num_threads();
        for(int i = 0; i < 300; i++) {
            cout << "-";
        }
        cout << endl;
        //cout << temp << ", num threads = " << omp_get_num_threads() << endl;

        if(root->left != nullptr) {
            preOrderTraversal(root->left, level + 1, sum);
        }

        if(root->right != nullptr) {
            preOrderTraversal(root->right, level + 1, sum);
        }
        
    //}
}




treeNode* insert(treeNode* root, double x) {
    if (root == nullptr) {
        treeNode* temp = new treeNode();
        temp->val = x;
        return temp;
    }

    if (x < root->val) {
        root->left = insert(root->left, x);
    } else if (x > root->val) {
        root->right = insert(root->right, x);
    }

    return root;
}
/*
treeNode* vecToTree(vector<double>& nums) {
    if(nums.size()==0)return NULL;
    if(nums.size()==1) return new treeNode(nums[0]);
    int middle = nums.size()/2;
    treeNode* root = new treeNode();
    root->val = nums[middle];
    vector<double> leftsub(nums.begin(), nums.begin() + middle);
    vector<double> rightsub(nums.begin() + middle+1, nums.end());
    root->left = vecToTree(leftsub);
    root->right = vecToTree(rightsub);
    return root;
}*/

void breadthFirstSearch(treeNode* root) {
    if (root == nullptr) {
        return;
    }

    queue<treeNode*> q;
    q.push(root);

    while (!q.empty()) {
        treeNode* current = q.front();
        q.pop();

        std::cout << current->val << " ";

        if (current->left != nullptr) {
            q.push(current->left);
        }

        if (current->right != nullptr) {
            q.push(current->right);
        }
    }
}





int main(int argc, char **argv) {
    //Usage: recursion_tree.out {n_threads, p1, p2, n, win_w, win_h}
    treeNode* root = new treeNode();
    treeNode* root2 = new treeNode();
    root->val = 0;
    vector<double> values;
    double start;
    double end;
    srandom(time(NULL));
    omp_set_num_threads(24);
    /*int nThreads = stoi(argv[1]);
    double p1 = stod(argv[2]);
    double p2 = stod(argv[3]);
    int n = stoi(argv[4]);
    int width = stoi(argv[5]);
    int height = stoi(argv[6]);
    cout << nThreads << ", " << p1 << ", " << p2 << endl;
    cout << "width = " << width << ", height = " << height << endl; */
    

    double min = 0;
    double max = 1;
    const long max_rand = 100000L;
    double sum = 0;


    for(int i = 0; i < 100000; i++) {
        values.push_back(min + (max - min) * (random() % max_rand) / max_rand);
        //cout << values[i] << endl;
    }
    
    root->val = values[0];
    for(int i = 1; i < values.size(); i++) {
       insert(root, values[i]);
       insert(root2, values[i]);
    }
    
    cout << "In-order traversal of the binary tree: " << endl;
    cout << "Starting num threads = " << omp_get_num_threads << endl;
    start = omp_get_wtime();
    omp_set_num_threads(12);
    #pragma omp parallel //default(none) shared(sum)
    {
        omp_set_num_threads(24);
        #pragma omp single 
        {
            cout << "Starting num threads = " << omp_get_num_threads() << endl;
            #pragma omp task //default(none) shared(sum)
            {
                sum = sum_recursive(root, 0);
                //preOrderTraversalOMP(root, 0, sum);
            }
        }
        /*#pragma omp task
        {
            preOrderTraversalOMP(root->right, 0);
        }*/
    }
    end = omp_get_wtime();
    cout << endl << "Parallel time = " << end - start << ", sum = " << sum << endl;

    start = omp_get_wtime();
    sum = sum_recursive(root, 0);
    //preOrderTraversal(root2, 0, sum);
    end = omp_get_wtime();
    cout << endl << "Normal time = " << end - start << ", sum = " << sum << endl;

    //breadthFirstSearch(root);
    //cout << "test print:" << endl;
    //printTree(root);
    cout << endl;


    return 0;
}