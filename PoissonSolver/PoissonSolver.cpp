#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
using namespace std;
using Func = function<double(double,double)>;
const double PI = 3.1415926535897932384626;
Func f = [](double x,double y){
        return 2*PI*PI*sin(PI*x)*sin(PI*y);
};
Func g = [](double x,double y){
        return 0.0;
};
Func exact_solution = [](double x,double y){
        return sin(PI*x)*sin(PI*y);
};
class PoissonSolver{
    private:
        int N;//剖分数，格点数=N+1
        double h;
        Func f;
        Func g;
        vector<vector<double>> u;
        
        double dotProduct (const vector<double>& a, const vector<double>& b){
            double result = 0.0;
            for(int i = 0;i < a.size();i++){
                result += a[i]*b[i];
            }
            return result;
        }
        void matrixVecProduct(const vector<double>& x, vector<double>& Ax){
            int n = (N-1)^2;
            fill(Ax.begin(), Ax.end(), 0.0);
            for (int i = 1;i < N;i++){
                for (int j = 1;j < N;j++){
                    int idx = (i-1)*(N-1)+j-1;
                    Ax[idx] = 4.0*x[idx];
                    if(i > 1){
                        int up_idx = (i-2)*(N-1)+j-1;
                        Ax[idx] -= x[up_idx];
                    }
                    if(i < N-1){
                        int down_idx = i*(N-1)+j-1;
                        Ax[idx] -= x[down_idx];
                    }
                    if(j > 1){
                        int left_idx = (i-1)*(N-1)+j-2;
                        Ax[idx] -= x[left_idx];
                    }
                    if(j < N-1){
                        int right_idx = (i-1)*(N-1)+j;
                        Ax[idx] -= x[right_idx];
                    }
                }
            }
        }
        void CG(const vector<double>& b, vector<double>& x, double tol, int max_iter){
            int n = b.size();
            vector<double> r = b;
            vector<double> p = r;
            vector<double> Ap(n);
            double b_norm = sqrt(dotProduct(b,b));
            if (b_norm < 1e-10) b_norm = 1.0;
            double rold = dotProduct(r, r);
            
            for (int k = 0; k < max_iter;k++){
                matrixVecProduct(p,Ap);
                double pAp = dotProduct(p,Ap);
                double alpha = rold/pAp;
                for (int i=0;i<n;i++){
                    x[i] +=alpha*p[i];
                    r[i] -=alpha*Ap[i];
                }
                double rnew = dotProduct(r,r);
                double r_norm = sqrt(rnew);
                if(r_norm/b_norm < tol) break;
                double beta = rnew/rold;
                for (int i=0;i<n;i++){
                    p[i] = r[i] + beta*p[i];
                }
                rold = rnew;
            }
        }
    public:
        
        PoissonSolver(int grid_size, Func source, Func boundary):N(grid_size),f(source),g(boundary){
            h = 1.0/N;
            u.resize(N+1,vector<double>(N+1,0.0));
            applyBC();
        }
        void applyBC(){
            for(int i=0;i<=N;i++){
                for(int j=0;j<=N;j++){
                    if(i==0 || i==N || j==0 || j==N){
                        double x = i*h;
                        double y = j*h;
                        u[i][j] = g(x,y);
                    }
                }
            }
        }
        void RHS(vector<double>& b){
            int n = (N-1)*(N-1);
            b.assign(n,0.0);
            for(int i=1;i<N;i++){
                for(int j=1;j<N;j++){
                    int idx = (i-1)*(N-1)+j-1;
                    double x=i*h;
                    double y=j*h;
                    b[idx] = h*h*f(x,y);
                    if(i==1) b[idx] += u[0][j];
                    if(i==N-1) b[idx] += u[N][j];
                    if(j==1) b[idx] += u[i][0];
                    if(j==N-1) b[idx] += u[i][N];
                }
            }
        }
        void solve(){
            int n = (N-1)*(N-1);
            if(n==0) return;
            vector<double> b;
            RHS(b);
            vector<double> x(n,0.0);
            CG(b,x,1e-8,1000000);
            for(int i=1;i<N;i++){
                for(int j=1;j<N;j++){
                    int idx = (i-1)*(N-1)+j-1;
                    u[i][j] = x[idx];
                }
            }
        }
        const vector<vector<double>>& getSolution() const {
            return u;
        }
};
int main(){
    int N = 10000;
    double h = 1.0/N;
    PoissonSolver solver(N,f,g);
    solver.solve();
    auto u = solver.getSolution();
    vector<double> x_grid(N+1),y_grid(N+1);
    for(int i=0;i<=N;i++){
        x_grid[i] = i*h;
        y_grid[i] = i*h;
    }
    vector<vector<double>> exact_u,error;
    exact_u.resize(N+1,vector<double>(N+1,0.0));
    error.resize(N+1,vector<double>(N+1,0.0));
    double max_error = 0.0;
    for (int i=0;i<=N;i++){
        for (int j=0;j<=N;j++){
            exact_u[i][j] = exact_solution(x_grid[i], y_grid[j]);
            error[i][j] = abs(u[i][j] - exact_u[i][j]);
            max_error = max(max_error, error[i][j]);
        }
    }
    cout << "Max error: " << max_error << endl;
    system("pause");
    return 0;
}
