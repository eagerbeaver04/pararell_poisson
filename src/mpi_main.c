#include "utils/utils.h"
#include <assert.h>
#include <mpi.h>
#include <stdio.h>

Matrix generate_five_diag(size_t xn, size_t yn)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t n, nn;
    if(mul_overflow_size_t(xn, yn, &n) || mul_overflow_size_t(n, n, &nn))
    {
        fprintf(stderr, "Слишком большая сетка (переполнение).\n");
        exit(EXIT_FAILURE);
    }

    Matrix A;

    // Only root process creates the actual matrix
    if(rank == 0)
    {
        A = create_matrix(n, n);
        printf("Main process: Created matrix %zux%zu\n", n, n);
    }

    // Broadcast dimensions to all processes
    MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xn, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Calculate work distribution
    int rows_per_process = n / size;
    int remainder = n % size;

    int start_row =
        rank * rows_per_process + (rank < remainder ? rank : remainder);
    int end_row = start_row + rows_per_process + (rank < remainder ? 1 : 0);
    int my_row_count = end_row - start_row;

    if(rank == 0)
    {
        // Main process computes its rows directly
        printf("Main process: Computing %d rows\n", my_row_count);
        for(size_t i = start_row; i < end_row; ++i)
        {
            A.data[i][i] = -4.0;
            if(i >= xn)
                A.data[i][i - xn] = 1.0;
            if(i + xn < n)
                A.data[i][i + xn] = 1.0;
            if((i % xn) != 0)
                A.data[i][i - 1] = 1.0;
            if((i % xn) != xn - 1)
                A.data[i][i + 1] = 1.0;
        }

        // Receive sparse data from workers
        for(int worker = 1; worker < size; worker++)
        {
            int worker_start = worker * rows_per_process +
                               (worker < remainder ? worker : remainder);
            int worker_end =
                worker_start + rows_per_process + (worker < remainder ? 1 : 0);
            int worker_row_count = worker_end - worker_start;

            if(worker_row_count > 0)
            {
                // Each worker sends: worker_row_count * 10 elements
                // Format: [col1, col2, col3, col4, col5, val1, val2, val3,
                // val4, val5] for each row
                double* sparse_data =
                    (double*)malloc(worker_row_count * 10 * sizeof(double));

                MPI_Recv(sparse_data, worker_row_count * 10, MPI_DOUBLE, worker,
                         0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Reconstruct matrix from sparse data
                for(int i = 0; i < worker_row_count; i++)
                {
                    int row = worker_start + i;
                    double* row_data = &sparse_data[i * 10];

                    // First 5 elements: column positions
                    // Next 5 elements: values
                    for(int j = 0; j < 5; j++)
                    {
                        int col = (int)row_data[j];     // Column index
                        double value = row_data[j + 5]; // Value

                        if(col >= 0)
                        { // Valid column (not -1)
                            A.data[row][col] = value;
                        }
                    }
                }

                free(sparse_data);
                printf("Main process: Received %d rows from worker %d\n",
                       worker_row_count, worker);
            }
        }

        printf("Main process: Matrix assembly complete\n");
    }
    else
    {
        // Workers compute and send sparse representation
        if(my_row_count > 0)
        {
            printf("Worker %d: Computing %d rows in sparse format\n", rank,
                   my_row_count);

            // Allocate buffer: 10 elements per row (5 positions + 5 values)
            double* sparse_data =
                (double*)malloc(my_row_count * 10 * sizeof(double));

            for(int i = 0; i < my_row_count; i++)
            {
                size_t row = start_row + i;
                double* row_data = &sparse_data[i * 10];

                // Initialize all positions to -1 and values to 0
                for(int j = 0; j < 10; j++)
                {
                    row_data[j] =
                        (j < 5) ? -1.0 : 0.0; // Positions = -1, Values = 0
                }

                // Center (always exists)
                row_data[0] = (double)row; // Position
                row_data[5] = -4.0;        // Value

                // Up (i - xn)
                if(row >= xn)
                {
                    row_data[1] = (double)(row - xn); // Position
                    row_data[6] = 1.0;                // Value
                }

                // Down (i + xn)
                if(row + xn < n)
                {
                    row_data[2] = (double)(row + xn); // Position
                    row_data[7] = 1.0;                // Value
                }

                // Left (i - 1)
                if((row % xn) != 0)
                {
                    row_data[3] = (double)(row - 1); // Position
                    row_data[8] = 1.0;               // Value
                }

                // Right (i + 1)
                if((row % xn) != xn - 1)
                {
                    row_data[4] = (double)(row + 1); // Position
                    row_data[9] = 1.0;               // Value
                }
            }

            // Send sparse data to main process
            MPI_Send(sparse_data, my_row_count * 10, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
            free(sparse_data);

            printf("Worker %d: Sent sparse data for %d rows to main process\n",
                   rank, my_row_count);
        }
        else
        {
            printf("Worker %d: No rows to compute\n", rank);
        }
    }

    return A; // Only meaningful for rank 0
}

Matrix cholesky(Matrix* A, int n)
{
    Matrix L = create_matrix(n, n);
    assert(A->cols == A->rows);
    assert(A->cols == n);

    for(int j = 0; j < n; j++)
    {
        double s = 0;
        for(int k = 0; k < j; k++)
        {
            s += L.data[j][k] * L.data[j][k];
        }
        L.data[j][j] = sqrt(A->data[j][j] - s);
#pragma omp parallel for
        for(int i = j + 1; i < n; i++)
        {
            double s = 0;
            for(int k = 0; k < j; k++)
            {
                s += L.data[i][k] * L.data[j][k];
            }
            L.data[i][j] = (1.0 / L.data[j][j] * (A->data[i][j] - s));
        }
    }
    return L;
}

Vector solve_gauss_reverse(Matrix* U, Vector* b)
{
    int n = U->rows;
    Vector x = create_vector(n);

    for(int i = n - 1; i >= 0; i--)
    {
        double sum = b->data[i];

        for(int j = i + 1; j < n; j++)
        {
            sum -= U->data[i][j] * x.data[j];
        }

        if(fabs(U->data[i][i]) < 1e-12)
        {
            printf("Ошибка: нулевой диагональный элемент в строке %d!\n", i);
            free_vector(x);
            exit(EXIT_FAILURE);
        }

        x.data[i] = sum / U->data[i][i];
    }

    return x;
}

Vector solve_gauss_forward(Matrix* L, Vector* b)
{
    int n = L->rows;
    Vector x = create_vector(n);

    for(int i = 0; i < n; i++)
    {
        double sum = b->data[i];

        for(int j = 0; j < i; j++)
        {
            sum -= L->data[i][j] * x.data[j];
        }

        if(fabs(L->data[i][i]) < 1e-12)
        {
            printf("Ошибка: нулевой диагональный элемент в строке %d!\n", i);
            free_vector(x);
            exit(EXIT_FAILURE);
        }

        x.data[i] = sum / L->data[i][i];
    }

    return x;
}

Vector solve_gauss(Matrix* L, Matrix* U, Vector* b)
{
    if(L->rows != U->rows || L->rows != b->size)
    {
        printf("Ошибка: несовместимые размеры в solve_gauss!\n");
        exit(EXIT_FAILURE);
    }

    Vector y = solve_gauss_forward(L, b);

    Vector x = solve_gauss_reverse(U, &y);

    free_vector(y);

    return x;
}

Vector pcgPreconditioned(Matrix* A, Vector* b, Vector* xs, double err,
                         double* relres, int* iter, Matrix* P1, Matrix* P2)
{
    int k = 0;
    size_t n = A->cols;
    Vector x = copy_vector(xs);
    Vector r = residue(b, A, &x);

    Vector z = solve_gauss(P1, P2, &r);

    Vector p = copy_vector(&z);

    double r0norm = second_norm(&r);

    Vector current_ = create_vector(n);
    Vector newR = create_vector(n);
    Vector newZ = create_vector(n);
    Vector q = create_vector(n);
    while(second_norm(&r) / r0norm > err && k < 1000)
    {
        ++k;

        free_vector(q);
        q = matrix_vector_mult(A, &p);

        double pq = dot_product(&p, &q);

        double a = aKbyPQ(&z, &q, pq);

        free_vector(current_);
        current_ = scalar_vector_mult(&p, a);

        add_vector_self(&x, &current_);

        free_vector(current_);
        current_ = scalar_vector_mult(&q, a);

        free_vector(newR);
        newR = sub_vector(&r, &current_);

        free_vector(newZ);
        newZ = solve_gauss(P1, P2, &newR);

        double b = dot_product(&newZ, &newR) / dot_product(&z, &r);

        free_vector(z);
        z = newZ;

        free_vector(r);
        r = newR;

        free_vector(current_);
        current_ = scalar_vector_mult(&p, b);

        free_vector(p);
        p = add_vector(&r, &current_);
    }
    *iter = k;

    free_vector(z);
    free_vector(r);
    free_vector(p);
    free_vector(current_);
    free_vector(newR);
    free_vector(newZ);
    free_vector(q);
    return x;
}

void errByEpsPcgChol(double a, double b, double c, double d, double h)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Vector x = linspace(a, b, (b - a) / h + 1);
    Vector y = linspace(c, d, (d - c) / h + 1);
    Matrix A = generate_five_diag(x.size - 2, y.size - 2);

    if(rank == 0)
    {
        Matrix L = cholesky(&A, A.cols);
        Vector us = uForXY(&x, &y);
        Vector B = F(&x, &y);
        int n = (x.size - 2) * (y.size - 2);

        printf("%zu, %zu, %i, %i, %i", A.cols, A.rows, B.size, x.size, y.size);
        scalar_mul_self(A, -1);

        scalar_vector_mult_self(&B, -1);

        printf("----\n");

        FILE* file = fopen("pcgCholErr.txt", "w");

        Matrix Lt = transpose(L);
        Vector zeros = create_vector(n);
        for(int i = 1; i <= 10; ++i)
        {
            double eps = pow(10, -i);
            double relres = 0;
            int count = 0;

            Vector Sol = pcgPreconditioned(&A, &B, &zeros, eps, &relres, &count,
                                           &L, &Lt);
            printf(", iter = %i\n", count);
            double max = vectors_max_diff(&us, &Sol);
            fprintf(file, "%.15f %.15f\n", eps, max);
            printf("%.15f %.15f\n", eps, max);
            if(i == 100)
            {
                for(int j = 0; j < us.size; ++j)
                {
                    printf("%.5f %.5f\n", us.data[j], Sol.data[j]);
                }
            }
            free_vector(Sol);
        }
        fclose(file);

        free_vector(us);
        free_vector(B);

        free_matrix(&Lt);
        free_vector(zeros);
        free_matrix(&L);
    }
    free_vector(x);
    free_vector(y);
    if(rank == 0)
    {
        free_matrix(&A);
    }
}

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0)
    {
        printf("Starting MPI program with %d processes\n", size);
    }

    double a = 0;
    double b = 1.625;
    double c = 0;
    double d = 1.625;
    double h = 0.025;
    printf("----------------------------------------------------\n");

    errByEpsPcgChol(a, b, c, d, h);
    printf("---\n");

    MPI_Finalize();

    return 0;
}