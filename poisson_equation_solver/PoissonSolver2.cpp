#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
// 周期边界条件下求解 -\Delta u + u = f
using namespace std;
using Func = function<double(double,double)>;
const double PI = 3.1415926535897932384626;
Func f = [](double x,double y){
        return (8*PI*PI+1)  *sin(2*PI*x)*sin(2*PI*y);
};
Func g = [](double x,double y){
        return 0.0;
};
Func exact_solution = [](double x,double y){
        return sin(2*PI*x)*sin(2*PI*y);
};
class PoissonSolver{
    private:
        int N;//剖分数，格点数=N+1
        double h=1.0/N;
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
            int n = N^2;
            double h = 1.0/N;
            fill(Ax.begin(), Ax.end(), 0.0);
            for (int i = 0;i < N;i++){
                for (int j = 0;j < N;j++){
                    int idx = i*N+j;
                    Ax[idx] = (4.0+h*h)*x[idx];
                    if (i == 0 ) Ax[idx] = Ax[idx] - x[idx+N] - x[idx+N*(N-1)];
                    else if (i == N-1) Ax[idx] = Ax[idx] - x[idx-N*(N-1)] - x[idx-N];
                    else Ax[idx] = Ax[idx] - x[idx-N] - x[idx+N];
                    if (j == 0) Ax[idx] = Ax[idx] - x[idx+1] - x[idx+N-1];
                    else if (j == N-1) Ax[idx] = Ax[idx] - x[idx-(N-1)] - x[idx-1];
                    else Ax[idx] = Ax[idx] - x[idx-1] - x[idx+1];   
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
            int n = N*N;
            b.assign(n,0.0);
            for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                    int idx = i*N+j;
                    double x=i*h;
                    double y=j*h;
                    b[idx] = h*h*f(x,y);
                    if(i==0) b[idx] += u[0][j];
                    if(j==0) b[idx] += u[i][0];
                }
            }
        }
        void solve(){
            int n = N*N;
            if(n==0) return;
            vector<double> b;
            RHS(b);
            vector<double> x(n,0.0);
            CG(b,x,1e-8,1000000);
            for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                    int idx = i*N+j;
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
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            exact_u[i][j] = exact_solution(x_grid[i], y_grid[j]);
            error[i][j] = abs(u[i][j] - exact_u[i][j]);
            max_error = max(max_error, error[i][j]);
        }
    }
    cout << "Max error: " << max_error << endl;
    system("pause");
    return 0;
}
