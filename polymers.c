/* Pivot Algorithm for Polymers
 * Author: Florian Seidler, Jonathan Zopes
 * Date: 14.02.2013
 * Issue: Computational Physics WT 12/13 , Bonn University */

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#define THREE_DIMENSIONS //Specify Dimension
#define MERSENNE_MAX 4294967295
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

using namespace std;

#ifdef TWO_DIMENSIONS
const int dimension=2;
#endif
#ifdef THREE_DIMENSIONS
const int dimension=3;
#endif
#ifdef FOUR_DIMENSIONS
const int dimension=4;
#endif

/*Tree nodes*/
struct saw_node {
    int n; //Number of Sites
    saw_node *left; //Left node pointer
    saw_node *right; //Right node pointer
    saw_node *parent; //Parent node pointer
    int** q; //Transformation Matrix
    int X[dimension]; //End Vector
    int B[2*dimension]; //Bounding box
};

/*-------------------*/
/* Random Number Gen */
/*-------------------*/

/* initializing the array with a NONZERO seed */
void
sgenrand(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

unsigned long
genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y;
}

/* Returns an integer in the interval [a,b) selected uniformly at random */
int random_integer_uniform(int a, int b) {
    if(b<=a) {
        cout<<"b is smaller than a, error!"<<endl;
        return 0;
    }
    double r = b - a + 1;
    return a + (int)(r * genrand()/MERSENNE_MAX);
}

int random_symmetry() {
    int sym;
    #ifdef TWO_DIMENSIONS
    int sym_num = 7;
    #endif // TWO_DIMENSIONS
    #ifdef THREE_DIMENSIONS
    int sym_num = 47;
    #endif // THREE_DIMENSIONS
    #ifdef FOUR_DIMENSIONS
    int sym_num = 383;
    #endif // FOUR_DIMENSIONS
    sym=random_integer_uniform(2,sym_num);
    return sym;
}

/*-------------------*/
/*Internal Operations*/
/*-------------------*/

/*Creates new quadratic matrix*/
int** new_sqmat() {
    int** q;
    q = new int*[dimension];
    for(int i = 0; i < dimension; i++) {
        q[i] = new int[dimension];
    }
    return q;
}

/*Deletes quadratic matrix*/
void delete_sqmat(int** q) {
    for(int i = 0; i < dimension; i++) delete[] q[i];
    delete[] q;
}

/*Sets single array to constant value*/
void set_single_array(int c, int* array) {
    for(int i=0;i<dimension;i++) {
        array[i]=c;
    }
}

/*Sets double array to constant value*/
void set_array(int c,int** array) {
    for(int i=0;i<dimension;i++) {
        for(int j=0;j<dimension;j++) {
            array[i][j]=c;
        }
    }
}

/*Copy array to array elementwise*/
void copy_array(int** input, int** result) {
    for(int i=0; i<dimension; i++)
    {
        for(int j=0; j<dimension; j++)
        {
            result[i][j]=input[i][j];
        }
    }
}

/*Copy single array to single array elementwise*/
void copy_single_array(int* input, int* result) {
    for(int i=0; i<dimension; i++) {
        result[i]=input[i];
    }
}

/*Set bounding box elements to constant value*/
void set_box(int c, int* array) {
    for(int i=0;i<2*dimension;i++) {
        array[i]=c;
    }
}

/*Matrix multiplication of two quadratic matrices*/
void mat_mul(int** left_op, int** right_op, int** product) {
    for(int i = 0; i<dimension; i++) {
        for(int j = 0; j<dimension; j++) {

            product[i][j] = 0; //product does not need to be initialized to 0.

            for(int k = 0; k<dimension; k++) {
                product[i][j] += left_op[i][k]*right_op[k][j];
            }
        }
    }
}

/*Matrix with vector multiplication*/
void mat_vec_mul(int** mat, int* vec, int* result) {
    for(int i=0; i<dimension; i++) {

        result[i] = 0; //result does not need to be initialized to 0.

        for(int j=0; j<dimension; j++) {
            result[i]+=mat[i][j]*vec[j];
        }
    }
}

/*Transpose Matrix method*/
void mat_transpose(int** mat, int**result) {

    for(int i=0; i<dimension; i++) {
        for(int j=0; j<dimension; j++) {
            result[i][j]=mat[j][i];
        }
    }

}

/*Sets up the transformation matrix specified by the number q.*/
void create_transformation_matrix(int q,int** result) {

    set_array(0,result); //initialize to zero

    #ifdef TWO_DIMENSIONS
    //in total 8 symmetry operations
    switch(q) {
        case 0: result[0][0]++; result[1][1]++; break; //identity -> do nothing
        case 1: result[0][1]++; result[1][0]++;  break; //reflect in 45deg
        case 2: result[0][1]--; result[1][0]--; break; //reflect in -45deg
        case 3: result[0][1]++; result[1][0]--; break; //rotate by 90deg clockwise
        case 4: result[0][0]--; result[1][1]--; break; //rotate by 180deg clockwise
        case 5: result[0][1]--; result[1][0]++; break; //rotate by 270deg clockwise
        case 6: result[0][0]++; result[1][1]--; break; //reflect at x-axis
        case 7: result[0][0]--; result[1][1]++; break; //reflect at y-axis
        default:
        cout<<"Wrong symmetry operation: "<<q<<endl;
        break;
    }
    #endif

    #ifdef THREE_DIMENSIONS
    //in total 48 symmetry operations
    switch(q%6) {
        case 0:   result[0][0]++; result[1][1]++; result[2][2]++;  break; // identity
        case 1:   result[0][1]++; result[1][0]++; result[2][2]++;  break; // (xy)
        case 2:   result[0][0]++; result[1][2]++; result[2][1]++;  break; // (yz)
        case 3:   result[0][2]++; result[1][1]++; result[2][0]++;  break; // (xz)
        case 4:   result[0][2]++; result[1][0]++; result[2][1]++;  break; // (xyz)
        case 5:   result[0][1]++; result[1][2]++; result[2][0]++;  break; // (xzy)
        default:
        cout<<"Wrong symmetry operation: "<<q<<endl;
        break;
    }
    switch(q/6) {
        case 0:   break; // +++
        case 1:   result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; break; // ++-
        case 2:   result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; break; // +-+
        case 3:   result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1;
                  result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; break; // +--
        case 4:   result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; break; // -++
        case 5:   result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1;
                  result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; break; // -+-
        case 6:   result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1;
                  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; break; // --+
        case 7:   result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1;
                  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1;
                  result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; break; // ---
        default:
        cout<<"Wrong symmetry operation: "<<q<<endl;
        break;
    }
    #endif

    #ifdef FOUR_DIMENSIONS
    switch(q%24) {
        case 0 :   result[0][0]++;   result[1][1]++;   result[2][2]++;   result[3][3]++;  break;
        case 1 :   result[0][0]++;   result[1][1]++;   result[2][3]++;   result[3][2]++;  break;
        case 2 :   result[0][0]++;   result[1][2]++;   result[2][1]++;   result[3][3]++;  break;
        case 3 :   result[0][0]++;   result[1][2]++;   result[2][3]++;   result[3][1]++;  break;
        case 4 :   result[0][0]++;   result[1][3]++;   result[2][1]++;   result[3][2]++;  break;
        case 5 :   result[0][0]++;   result[1][3]++;   result[2][2]++;   result[3][1]++;  break;
        case 6 :   result[0][1]++;   result[1][0]++;   result[2][3]++;   result[3][2]++;  break;
        case 7 :   result[0][1]++;   result[1][0]++;   result[2][2]++;   result[3][3]++;  break;
        case 8 :   result[0][1]++;   result[1][2]++;   result[2][0]++;   result[3][3]++;  break;
        case 9 :   result[0][1]++;   result[1][2]++;   result[2][3]++;   result[3][0]++;  break;
        case 10:   result[0][1]++;   result[1][3]++;   result[2][0]++;   result[3][2]++;  break;
        case 11:   result[0][1]++;   result[1][3]++;   result[2][2]++;   result[3][0]++;  break;
        case 12:   result[0][2]++;   result[1][0]++;   result[2][1]++;   result[3][3]++;  break;
        case 13:   result[0][2]++;   result[1][0]++;   result[2][3]++;   result[3][1]++;  break;
        case 14:   result[0][2]++;   result[1][1]++;   result[2][0]++;   result[3][3]++;  break;
        case 15:   result[0][2]++;   result[1][1]++;   result[2][3]++;   result[3][0]++;  break;
        case 16:   result[0][2]++;   result[1][3]++;   result[2][0]++;   result[3][1]++;  break;
        case 17:   result[0][2]++;   result[1][3]++;   result[2][1]++;   result[3][0]++;  break;
        case 18:   result[0][3]++;   result[1][0]++;   result[2][1]++;   result[3][2]++;  break;
        case 19:   result[0][3]++;   result[1][0]++;   result[2][2]++;   result[3][1]++;  break;
        case 20:   result[0][3]++;   result[1][1]++;   result[2][0]++;   result[3][2]++;  break;
        case 21:   result[0][3]++;   result[1][1]++;   result[2][2]++;   result[3][0]++;  break;
        case 22:   result[0][3]++;   result[1][2]++;   result[2][0]++;   result[3][1]++;  break;
        case 23:   result[0][3]++;   result[1][2]++;   result[2][1]++;   result[3][0]++;  break;
        default:
        cout<<"Wrong symmetry operation: "<<q<<endl;
        break;
    }
    switch(q/24) {
        case 0: break;
        case 1:  result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // +++-
        case 2:  result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1; break; // ++-+
        case 3:  result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // ++--
        case 4:  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1; break; // +-++
        case 5:  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // +-+-
        case 6:  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1; break; // +--+
        case 7:  result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // +---
        case 8:  result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1; break; // -+++
        case 9:  result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // -++-
        case 10: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1; break; // -+-+
        case 11: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // -+--
        case 12: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1; break; // --++
        case 13: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // --+-
        case 14: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1; break; // ---+
        case 15: result[0][0]*=-1; result[0][1]*=-1; result[0][2]*=-1; result[0][3]*=-1;
                 result[1][0]*=-1; result[1][1]*=-1; result[1][2]*=-1; result[1][3]*=-1;
                 result[2][0]*=-1; result[2][1]*=-1; result[2][2]*=-1; result[2][3]*=-1;
                 result[3][0]*=-1; result[3][1]*=-1; result[3][2]*=-1; result[3][3]*=-1; break; // ----
        default:
        cout<<"Wrong symmetry operation: "<<q<<endl;
        break;
   }
    #endif

}

/*Adds vectors together*/
void addX(int* point1, int* point2, int* result) {

    for(int i=0; i<dimension; i++) {
        result[i]=point1[i]+point2[i];
    }

}

/*Translation by vector X applied to bounding box B.*/
void transB(int* X, int* B, int* result) {

    for(int i=0; i<dimension; i++) {
        result[i]=B[i]+X[i];
        result[i+dimension]=B[i+dimension]+X[i];
    }

}

/* Fuses the boxes B1 and B2.*/
void fuseB(int* B1, int* B2, int* result) {

    for(int i=0;i<dimension;i++) {
        result[i]=fmin(B1[i],B2[i]);
        result[i+dimension]=fmax(B1[i+dimension],B2[i+dimension]);
    }

}

/* Checks bounding boxes for overlap. */
bool overlapB(int* B1, int* B2) {

    for(int i=0; i<dimension; i++) {
        if(fmax(B1[i],B2[i])>fmin(B1[i+dimension],B2[i+dimension])) return false;
    }

    return true;

}


void printB(int* B) {

    for(int i=0;i<dimension;i++) {
        cout<<B[i]<<"\t";
    }
    cout<<endl;

    for(int i=dimension;i<2*dimension;i++) {
        cout<<B[i]<<"\t";
    }
    cout<<endl;

}

void printQ(int** q) {

    for(int i=0;i<dimension;i++) {
        cout<<"(";
        for(int j=0;j<dimension;j++) {
            cout<<q[i][j]<<" ";
        }
        cout<<")"<<endl;
    }
    cout<<endl;

}

void printX(int* X) {

    for(int i=0; i<dimension; i++) {
        cout<<X[i]<<endl;
    }

}

/*-----------------------------------------*/
/*Tree Structure and Operations on the Tree*/
/*-----------------------------------------*/

/*Construct saw nodes with preset values*/
struct saw_node* construct_saw_node(saw_node* parentL, int nL, int** qL, int* XL, int* BL) {

    saw_node* construct=new saw_node;
    construct->parent=parentL;
    construct->n=nL;
    construct->left=NULL;
    construct->right=NULL;
    construct->q = new_sqmat();

    for(int i=0; i<dimension; i++) {
        for(int j=0; j<dimension; j++) {
            construct->q[i][j]=qL[i][j];
        }
    }

    for(int i=0; i<dimension; i++) {
        construct->X[i]=XL[i];
    }

    for(int i=0; i<2*dimension; i++) {
        construct->B[i]=BL[i];
    }
    return construct;

}

/*Construct trivial trees*/
int trivial_tree(saw_node* root, saw_node* leaf) {

    int n;
    int X[dimension];
    int B[2*dimension];
    set_single_array(0,X);
    set_box(0,B);
    int **q=new_sqmat();
    create_transformation_matrix(0,q);
    saw_node* curr=root;
    n=curr->n;

    if (n>2) {
        n/=2;
        X[0]=n;
        B[0]=1; // Test to initialize the correct bounding box.
        B[dimension]=X[0];
        curr->left = construct_saw_node(curr,n,q,X,B);
        curr->right = construct_saw_node(curr,n,q,X,B);
        trivial_tree(curr->left,leaf);
        trivial_tree(curr->right,leaf);
    } else {
        curr->left=leaf;
        curr->right=leaf;
    }

    delete_sqmat(q);

    return 0;
}

/*Generates saw tree with n nodes.*/
struct saw_node *generate_saw_tree(int n) {

    if(n<=0){
    cout<<"invalid saw-tree size:"<<n<<endl;
    return NULL;
    }

    int D,L;
    int X[dimension];
    int B[2*dimension];
    set_single_array(0,X);
    set_box(0,B);
    X[0]=n;
    B[0]=1; // Test to initialize the correct bounding box.
    B[dimension]=n;
    int **q=new_sqmat();
    create_transformation_matrix(0,q);
    saw_node *root = construct_saw_node(NULL,n,q,X,B);
    X[0]=1;
    B[0]=1; // Test to initialize the correct bounding box.
    B[dimension]=1;
    saw_node *leaf = construct_saw_node(NULL,1,q,X,B);
    D=ceil(log((double)n)/log(2.0));
    saw_node* curr = root;

    for(L=1;L<D-1;L++) {
        if(n-rint(ldexp(1,D-L))>= rint(ldexp(1,D-L-1))) {
            X[0] = rint(ldexp(1,D-L));
            B[0] = 1; // Test to initialize the correct bounding box.
            B[dimension]=X[0];
            curr->left = construct_saw_node(curr,X[0],q,X,B);
            X[0] = n-rint(ldexp(1,D-L));
            B[0] = 1; // Test to initialize the correct bounding box.
            B[dimension]=X[0];
            curr->right = construct_saw_node(curr,X[0],q,X,B);
            //Bounding box initial
            trivial_tree(curr->left, leaf);
            curr=curr->right;
        } else {
            X[0] = n-rint(ldexp(1,D-L-1));
            B[0] = 1; // Test to initialize the correct bounding box.
            B[dimension]=X[0];
            curr->left = construct_saw_node(curr,X[0],q,X,B);
            X[0] = rint(ldexp(1,D-L-1));
            B[0] = 1; // Test to initialize the correct bounding box.
            B[dimension]=X[0];
            curr->right = construct_saw_node(curr,X[0],q,X,B);
            trivial_tree(curr->right, leaf);
            curr=curr->left;
        }
        n=curr->n;
    }

    switch(n) {
        case 1:     curr=leaf;
        case 2:     curr->left = leaf;
                    curr->right = leaf;
                    break;
        case 3:     X[0] = 2;
                    B[0] = 1; // Test to initialize the correct bounding box.
                    B[dimension]=X[0];
                    curr->left = construct_saw_node(curr,X[0],q,X,B);
                    curr->right = leaf;
                    curr->left->left = leaf;
                    curr->left->right = leaf;
                    break;
        case 4:     trivial_tree(curr, leaf);
                    break;
        default:    cout << "Error in tree generation." << endl;
                    break;
    }

    delete_sqmat(q);
    return root;

}

/*Clear tree structure and free memory*/
saw_node* clear_internal(saw_node* node) {

    saw_node* leaf = NULL;
    if(node->right || node->left){
        leaf = clear_internal(node->left);
        clear_internal(node->right);
    } else {
        return node;
    }

    delete_sqmat(node->q);
    delete node;
    return leaf;

}

/*Clear the tree structure method*/
void clear_saw_tree(saw_node* node) {

    saw_node* leaf;
    leaf= clear_internal(node);
    delete_sqmat(leaf->q);
    delete leaf;

}

/*Merge tree wl and wr in tree node w*/
void merge_saw_tree(saw_node* wl, saw_node* wr, saw_node* w) {

    //Allocate buffer for transformations
    int bufferX1[dimension];
    int bufferX2[dimension];
    int bufferY1[dimension];
    int bufferY2[dimension];
    int bufferB1[2*dimension];
    int bufferB2[2*dimension];

    //Set pointers
    w->left=wl;
    w->right=wr;

    //n=n_l+n_r
    w->n=wl->n+wr->n;

    //X=X_l + q*X_r
    //copy_array(qL,w->q);
    mat_vec_mul(w->q,wr->X,bufferX1);
    addX(wl->X,bufferX1,w->X);

    //B=B_l u (X_l + q*B_r)
    //B_r stored in 2 vectors which will be transformed
    for(int i=0; i<dimension;i++) {
        bufferX1[i]=w->right->B[i];
        bufferX2[i]=w->right->B[i+dimension];
    }
    set_single_array(0,bufferY1);
    set_single_array(0,bufferY2);
    mat_vec_mul(w->q,bufferX1,bufferY1);
    mat_vec_mul(w->q,bufferX2,bufferY2);

    for(int i=0; i<dimension;i++) {
        bufferB1[i]=fmin(bufferY1[i],bufferY2[i]);
        bufferB1[i+dimension]=fmax(bufferY1[i],bufferY2[i]);
    }

    //Translate bufferB1 by X_l and store in bufferB2
    transB(wl->X,bufferB1,bufferB2);

    //Fuse bufferB2 with B_l and store in B
    fuseB(bufferB2,wl->B,w->B);


}

/*Recursive print call*/
int rec_print(saw_node* node, saw_node* parent, bool right_side, saw_node** node_Arr, int node_Arr_1stfree, fstream& outfile) {

    int XL[dimension], XLL[dimension];

    if(right_side) {
        node_Arr[node_Arr_1stfree] = parent;
        node_Arr_1stfree++;
    }
    if(node->right || node->left) {//not leaf
        rec_print(node->left, node, false, node_Arr, node_Arr_1stfree, outfile);
        rec_print(node->right, node, true, node_Arr, node_Arr_1stfree, outfile);
    } else {
        copy(node->X, node->X+dimension, XL);
        for(int i = node_Arr_1stfree-1; i >= 0; i--) {
            mat_vec_mul(node_Arr[i]->q, XL, XLL);
            addX(XLL, node_Arr[i]->left->X, XL);
        }
        for(int k = 0; k < dimension; k++) outfile << XL[k] << " ";
        outfile << endl;
    }
    if(right_side) {
        node_Arr[node_Arr_1stfree-1] = NULL;
        node_Arr_1stfree--;
    }

    return 0;

}

/*Actual print function*/
int print_tree(saw_node* root, char* outfile_name) {

    saw_node** node_Arr = new saw_node*[root->n];

    for(int i = 0; i < root->n; i++) node_Arr[i] = NULL;
    fstream outfile(outfile_name, fstream::out);
    rec_print(root, root->parent, false, node_Arr, 0, outfile);
    outfile.close();
    delete[] node_Arr;
    return 0;

}

/*write an Array with spaces terminated by a tab*/
int write_Array(int Arr[], int length, ofstream& outfile) {

    for(int i = 0; i < length; i++){
        outfile << Arr[i] << " ";
    }
    outfile << "\t";

    return 0;
}

/*Save tree by writing node to a file recursively*/
int save_tree_rec(saw_node* node, ofstream& outfile) {

    //write data
    outfile << node->n << "\t";
    for(int i = 0; i < dimension; i++){
        write_Array(node->q[i], dimension, outfile);
    }
    write_Array(node->X, dimension, outfile);
    write_Array(node->B, 2*dimension, outfile);
    outfile << endl;

    //recursion
    if(node->right || node->left){//not leaf
        save_tree_rec(node->left, outfile);
        save_tree_rec(node->right, outfile);
    }

    return 0;
}

/*Actual save function*/
int save_tree(saw_node* root, char* outfile_name){

    ofstream outfile(outfile_name, ios::out);
    save_tree_rec(root, outfile);
    outfile.close();

    return 0;
}

/*read an Array seperated by spaces*/
int read_Array(int Arr[], int length, ifstream& infile) {

    for(int i = 0; i < length; i++){
        infile >> Arr[i];
    }

    return 0;
}

/*Load a tree recursively by reading a file*/
int load_tree_rec(saw_node* &node, saw_node* parentL, saw_node** leaf, ifstream& infile) {

    //read data into node
    node = new saw_node;
    node->parent = parentL;
    infile >> node->n;
    node->q = new_sqmat();
    for(int i = 0; i < dimension; i++){
        read_Array(node->q[i], dimension, infile);
    }
    read_Array(node->X, dimension, infile);
    read_Array(node->B, 2*dimension, infile);

    if(node->n > 1){//not leaf

        //recursion
        load_tree_rec(node->left, node, leaf, infile);
        load_tree_rec(node->right, node, leaf, infile);
    }
    else if(!*leaf) {
            node->parent = NULL;
            node->left   = NULL;
            node->right  = NULL;
            leaf = &node;
    }
    else {
        delete_sqmat(node->q);
        delete node;
        node = *leaf;
    }

    return 0;
}

/*Actual load function*/
int load_tree(saw_node* &root, char* infile_name) {

    ifstream infile("test.txt", ios::in);
    saw_node** leaf = new saw_node*;
    *leaf = NULL;

    load_tree_rec(root, NULL, leaf, infile);
    delete leaf;
    infile.close();

    return 0;
}

/*Left tree-rotation applied to w*/
void LR(saw_node* w) {

    int** qt=new_sqmat();
    saw_node* wt;
    wt=w->right;
    w->right=wt->right;
    wt->right=wt->left;
    wt->left=w->left;
    w->left=wt;

    //set parents
    w->left->left->parent=w->left;
    w->right->parent=w;

    //change symmetry transformations
    copy_array(w->q,qt);
    mat_mul(qt,w->left->q,w->q);
    copy_array(qt,w->left->q);

    //merge(wll,wlr,wl->q)
    merge_saw_tree(w->left->left,w->left->right,w->left);

    //free memory
    delete_sqmat(qt);

}

/*Right tree-rotation applied to w*/
void RR(saw_node* w) {

    int** qbuffer=new_sqmat();
    int** qt=new_sqmat();
    saw_node* wt;
    wt=w->left;
    w->left=wt->left;
    wt->left=wt->right;
    wt->right=w->right;
    w->right=wt;

    //set parents
    w->right->right->parent=w->right;
    w->left->parent=w;

    //change symmetry transformations
    copy_array(w->q,qt);
    copy_array(w->right->q,w->q);
    mat_transpose(w->q,qbuffer);
    mat_mul(qbuffer,qt,w->right->q);
    merge_saw_tree(w->right->left,w->right->right, w->right);

    //free memory
    delete_sqmat(qbuffer);
    delete_sqmat(qt);

}

/*Shuffle up applied to tree w*/
void shuffle_up(int n0, saw_node* w) {

    if(n0<w->left->n) {
        shuffle_up(n0,w->left);
        RR(w);
    }
    else if(n0>w->left->n) {
        shuffle_up(n0-w->left->n,w->right);
        LR(w);
    }

    return;

}

/*Shuffle down applied to tree w*/
void shuffle_down(saw_node* w) {

    int nt;
    nt=(w->n+1)/2; //integer division replaces floor function
    if(nt<w->left->n) {
        RR(w);
        shuffle_down(w->right);
    }
    else if(nt>w->left->n) {
        LR(w);
        shuffle_down(w->left);
    }

    return;

}


/*Intersection test function*/
bool intersect(int* xl_abs, int** ql_abs, saw_node* wl, int* xr_abs, int** qr_abs, saw_node* wr)
{
    bool result;
    //Allocate temporary memory
    int bufferX1[dimension];
    int bufferX2[dimension];
    int bufferY1[dimension];
    int bufferY2[dimension];
    int bufferB1[2*dimension];
    int bufferB2[2*dimension];
    int bufferB3[2*dimension];
    int** qbuffer=new_sqmat();

    //B_l transformed and translated and then stored in bufferB2
    for(int i=0; i<dimension;i++) {
        bufferX1[i]=wl->B[i];
        bufferX2[i]=wl->B[i+dimension];
    }
    mat_vec_mul(ql_abs,bufferX1,bufferY1);
    mat_vec_mul(ql_abs,bufferX2,bufferY2);
    for(int i=0; i<dimension;i++) {
        bufferB1[i]=fmin(bufferY1[i],bufferY2[i]);
        bufferB1[i+dimension]=fmax(bufferY1[i],bufferY2[i]);
    }

    //Translate bufferB1 by X_l and store in bufferB2
    transB(xl_abs,bufferB1,bufferB2);

    //B_r transformed and translated and then stored in bufferB3
    for(int i=0; i<dimension;i++) {
        bufferX1[i]=wr->B[i];
        bufferX2[i]=wr->B[i+dimension];
    }
    mat_vec_mul(qr_abs,bufferX1,bufferY1);
    mat_vec_mul(qr_abs,bufferX2,bufferY2);
    for(int i=0; i<dimension;i++) {
        bufferB1[i]=fmin(bufferY1[i],bufferY2[i]);
        bufferB1[i+dimension]=fmax(bufferY1[i],bufferY2[i]);
    }

    //Translate bufferB1 by X_l and store in bufferB3
    transB(xr_abs,bufferB1,bufferB3);

    //Check for overlap
    if(!overlapB(bufferB2,bufferB3)) {
        delete_sqmat(qbuffer);
        return false;
    }
    if(wl->n<=2 && wr->n<=2) {
        delete_sqmat(qbuffer);
        return true;
    }

    if(wl->n>=wr->n) {
        //ql_abs*X_ll and store in bufferX1
        mat_vec_mul(ql_abs,wl->left->X,bufferX1);
        //xl_abs+bufferX1 and store in bufferX2
        addX(xl_abs,bufferX1,bufferX2);
        //ql_abs*wl->q and store in qbuffer
        mat_mul(ql_abs,wl->q,qbuffer);
        if(intersect(bufferX2,qbuffer,wl->right,xr_abs,qr_abs,wr)) {
            delete_sqmat(qbuffer);
            return true;
        }
        delete_sqmat(qbuffer);
        return intersect(xl_abs,ql_abs,wl->left,xr_abs,qr_abs,wr);
    }
    else {
        if(intersect(xl_abs,ql_abs,wl,xr_abs,qr_abs,wr->left)){
            delete_sqmat(qbuffer);
            return true;
        }
        //qr_abs*X_rl and store in bufferX1
        mat_vec_mul(qr_abs,wr->left->X,bufferX1);
        //xr_abs+bufferX1 and store in bufferX2
        addX(xr_abs,bufferX1,bufferX2);
        //qr_abs*qr and store in qbuffer
        mat_mul(qr_abs,wr->q,qbuffer);
        result=intersect(xl_abs,ql_abs,wl,bufferX2,qbuffer,wr->right);
        delete_sqmat(qbuffer);
        return result;
    }
}

/*Attempt pivot operation*/
bool attempt_pivot(saw_node* w, int nt, int** qt) {

    //allocate memory
    int** qbuffer1=new_sqmat();
    int** qbuffer2=new_sqmat();
    int X[dimension];
    bool intersection;

    shuffle_up(nt,w);

    //q<=q*qt
    mat_mul(w->q,qt,qbuffer1);
    copy_array(qbuffer1,w->q);

    //Store identity matrix in qbuffer2 and origin vector in X
    create_transformation_matrix(0,qbuffer1);
    set_single_array(0,X);

    //intersect(0,id,wl,wl->X,q,wr)
    intersection=intersect(X,qbuffer1,w->left,w->left->X,w->q,w->right);

    //reset the initial walk parameters if the transformation lead to intersection
    if(intersection){
        mat_transpose(qt,qbuffer1);
        mat_mul(w->q,qbuffer1,qbuffer2);
        copy_array(qbuffer2,w->q);
    }
    shuffle_down(w);
    delete_sqmat(qbuffer1);
    delete_sqmat(qbuffer2);
    return !intersection;
}

/*Apply dimerize operation to tree w*/
void dimerize(saw_node* w) {
    int** qt=new_sqmat();
    int sym;
    double attempts=20*pow(w->n,1.2);
    int numberofattempts=(int)attempts;
    cout<<"Dimerizing... Dimerize Attempts:"<<numberofattempts<<endl;
    for(int i=0; i<numberofattempts;i++) {
        sym=random_symmetry();
        create_transformation_matrix(sym,qt);
        attempt_pivot(w,random_integer_uniform(1,w->n-1),qt);
        if(i%30000==0) {
            cout.precision(1);
            cout<< fixed <<(double)i/(double)numberofattempts*100.0<<"%"<<endl;
        }
    }
    delete_sqmat(qt);
    cout<<"Dimerize done!"<<endl;
}

/*Calculate the square distance to the origin*/
double calc_square_distance(saw_node* w) {
    double distance=0;
    int buffer1[dimension];
    int buffer2[dimension];
    int endvector[dimension];
    copy_single_array(w->left->X,buffer2);
    mat_vec_mul(w->q,w->right->X,buffer1);
    addX(buffer2,buffer1,endvector);
    for(int i=0;i<dimension; i++) {
        distance+=endvector[i]*endvector[i];
    }
    return distance;
}

/*Gather statistics*/
double gather_statistics(saw_node* w, int samplesize) {
    bool success;
    double sum=0;
    int sym;
    int** qt=new_sqmat();
    int i=0;
    while(i<samplesize){
        sym=random_symmetry();
        create_transformation_matrix(sym,qt);
        success=attempt_pivot(w,random_integer_uniform(1,w->n-1),qt);
        if(success) {
           sum+=calc_square_distance(w);
           i++;
        }
    }
    delete_sqmat(qt);
    return sum/(double)samplesize;
}


int main() {
    saw_node* tree;
    tree=generate_saw_tree(50000);
    sgenrand(time(NULL));
    dimerize(tree);
    save_tree(tree,"dimerized-50k.txt");
    clear_saw_tree(tree);
    return 0;

}
