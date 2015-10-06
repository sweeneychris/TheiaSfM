#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void eigs_sym_F77(MatrixXd &M, VectorXd &init_resid, int k, int m, double &time_used, double &prec_err);
void eigs_gen_F77(MatrixXd &M, VectorXd &init_resid, int k, int m, double &time_used, double &prec_err);
void eigs_sym_Cpp(MatrixXd &M, VectorXd &init_resid, int k, int m, double &time_used, double &prec_err);
void eigs_gen_Cpp(MatrixXd &M, VectorXd &init_resid, int k, int m, double &time_used, double &prec_err);


void print_header(std::string title)
{
    const int width = 80;
    const char sep = ' ';

    std::cout << std::endl << std::string(width, '=') << std::endl;
    std::cout << std::string((width - title.length()) / 2, ' ') << title << std::endl;
    std::cout << std::string(width, '-') << std::endl;

    std::cout << std::left << std::setw(12) << std::setfill(sep) << "mat_size";
    std::cout << std::left << std::setw(11) << std::setfill(sep) << "dataset";
    std::cout << std::left << std::setw(15) << std::setfill(sep) << "F77_time/ms";
    std::cout << std::left << std::setw(13) << std::setfill(sep) << "F77_error";
    std::cout << std::left << std::setw(15) << std::setfill(sep) << "C++_time/ms";
    std::cout << std::left << std::setw(13) << std::setfill(sep) << "C++_error";
    std::cout << std::endl;

    std::cout << std::string(width, '-') << std::endl;
}

void print_row(int n, int dataset, double time_f77, double err_f77, double time_cpp, double err_cpp)
{
    const char sep = ' ';

    std::cout << std::left << std::setw(12) << std::setfill(sep) << n;
    std::cout << std::left << std::setw(11) << std::setfill(sep) << dataset;
    std::cout << std::left << std::setw(15) << std::setfill(sep) << time_f77;
    std::cout << std::left << std::setw(13) << std::setfill(sep) << err_f77;
    std::cout << std::left << std::setw(15) << std::setfill(sep) << time_cpp;
    std::cout << std::left << std::setw(13) << std::setfill(sep) << err_cpp;
    std::cout << std::endl;
}

void print_footer()
{
    const int width = 80;
    std::cout << std::string(width, '=') << std::endl << std::endl;
}

void run_eigs_sym(int n_experiment, int n_replicate, int n, int k, int m)
{
    double time_f77, time_cpp;
    double err_f77, err_cpp;

    for(int i = 0; i < n_experiment; i++)
    {
        MatrixXd A = MatrixXd::Random(n, n);
        MatrixXd M = A.transpose() + A;

        VectorXd init_resid = VectorXd::Random(M.cols());
        init_resid.array() -= 0.5;
        init_resid = M * init_resid;

        for(int j = 0; j < n_replicate; j++)
        {
            eigs_sym_F77(M, init_resid, k, m, time_f77, err_f77);
            eigs_sym_Cpp(M, init_resid, k, m, time_cpp, err_cpp);
            print_row(n, i + 1, time_f77, err_f77, time_cpp, err_cpp);
        }
    }
}

void run_eigs_gen(int n_experiment, int n_replicate, int n, int k, int m)
{
    double time_f77, time_cpp;
    double err_f77, err_cpp;

    for(int i = 0; i < n_experiment; i++)
    {
        MatrixXd A = MatrixXd::Random(n, n);

        VectorXd init_resid = VectorXd::Random(A.cols());
        init_resid.array() -= 0.5;
        init_resid = A * init_resid;

        for(int j = 0; j < n_replicate; j++)
        {
            eigs_gen_F77(A, init_resid, k, m, time_f77, err_f77);
            eigs_gen_Cpp(A, init_resid, k, m, time_cpp, err_cpp);
            print_row(n, i + 1, time_f77, err_f77, time_cpp, err_cpp);
        }
    }
}

int main()
{
    srand(123);
    int n_experiment = 5;
    int n_replicate = 5;

    print_header("eigs_sym");
    run_eigs_sym(n_experiment, n_replicate, 100, 10, 20);
    run_eigs_sym(n_experiment, n_replicate, 1000, 10, 30);
    print_footer();

    print_header("eigs_gen");
    run_eigs_gen(n_experiment, n_replicate, 100, 10, 20);
    run_eigs_gen(n_experiment, n_replicate, 1000, 10, 30);
    print_footer();

    return 0;
}
