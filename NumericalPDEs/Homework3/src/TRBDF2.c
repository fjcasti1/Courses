/*
*				TRBDF2.c
*	Copyright 2016 by Carl L. Gardner.  All rights reserved.
*
*	To compile: cc -o TRBDF2 TRBDF2.c -lm
*	Can use: cc -Ofast -o TRBDF2 TRBDF2.c -lm
*	To run: ./TRBDF2 > OUT
*	Or:     ./TRBDF2 < in-TRBDF2 > OUT &
*
*	Solves the nonlinear diffusion equation
*			du/dt = F = d/dx (D(u) du/dx)
*	using TRBDF2.
*
*	Originally written in 5 files:
*	cons.c  diffuse.c  parab.c  pdecs.h  pmain.c
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// #include <pdecs.h>

typedef void	*POINTER;	// pointer to an unknown data type

#define EPSILON_M 2.2204460492503131e-16
#define sq(x)	( (x)*(x) )
#define YES	(1)
#define NO	(0)
#define TRUE	(1)
#define FALSE	(0)
#define ERROR	(-1)

#define ABS(x)		( ((x) >= 0) ? (x) : -(x) )
#define max(a,b)	( ((a) > (b)) ? (a) : (b) )
#define min(a,b)	( ((a) < (b)) ? (a) : (b) )

#define INT		((unsigned) sizeof(int))
#define FLOAT		((unsigned) sizeof(float))
#define DOUBLE		((unsigned) sizeof(double))
#define CHAR		((unsigned) sizeof(char))

// library function prototypes; functions at end of file
POINTER Alloc(unsigned N_bytes);
POINTER alloc_vector(int N, unsigned element_size);
POINTER alloc_matrix(int N_rows, int N_columns,	unsigned element_size);
void compute_norm(double v[], int N, double *norm);

/*
*				pdecs.h
*
*	Contains the structures & definitions for semiconductor process
*	simulations, & defines the system of units.
*
*			Computational Units
*
*	Fundamental units
*		L_scale = 0.1 micron = 10^-5 cm
*		t_scale = 1000 sec
*		u_scale = 10^20 cm^-3
*
*	Derived units
*		D_scale = 10^-13 cm^2/sec = L_scale^2/t_scale
*
*	In other words, if D_phys = 2*10^-13 cm^2/sec in physical units,
*			D = D_phys/D_scale = 2
*	in computational units.
*/

// define computational units
#define u_scale		1.e20
#define D_scale		1.e-13
#define L_scale		0.1
#define t_scale		1000

typedef struct {
	int Nx;
	int modes;	// number of unknowns
	double dx, xmin, xmax;
} GRID;

typedef struct {
	int step, REDO_STEP, LAST_STEP;
	double dt, dt_min, dt_max, dt_euler, t, t_max;
	double FACTOR, MIN_FACTOR, MAX_FACTOR; // dt = FACTOR*dt_euler
	double norm_change, norm_u;
} TIMESTEP;

// function prototypes
// cons.c
void check_conservation(GRID *grid, double u[]);
double compute_total(GRID *grid, double u[]);
// diffuse.c
void compute_D(GRID *grid, double u[], double D[]);
void compute_deltaD(GRID *grid, double u[], double D[], double deltaD[]);
// parab.c
void TRBDF2(GRID *grid, TIMESTEP *timestep, double u[]);
void TR_step(double dt, GRID *grid, TIMESTEP *timestep, double u[],
	double u_mid[], double F[], double F_mid[], double D[],
	double deltaD[], double *J[], double b[]);
void BDF2_step(double dt, GRID *grid, TIMESTEP *timestep, double u[],
	double u_mid[], double u_new[], double F_new[], double D[],
	double deltaD[], double *J[], double b[]);
void compute_F(GRID *grid, double u[], double D[], double F[]);
void compute_deltaF(GRID *grid, double u[], double D[], double deltaD[],
	double *J[]);
void compute_TR_J(double dt, GRID *grid, double *J[]);
void compute_BDF2_J(double dt, GRID *grid, double *J[]);
void compute_TR_b(double dt, GRID *grid, double u[], double u_mid[],
	double F[], double F_mid[], double b[]);
void compute_BDF2_b(double dt, GRID *grid, double u[], double u_mid[],
	double u_new[], double F_new[], double b[]);
void update_states(GRID *grid, TIMESTEP *timestep, double u[], double b[]);
void compute_dt(GRID *grid, TIMESTEP *timestep, double F[], double F_mid[],
	double F_new[]);
void compute_first_dt(GRID *grid, TIMESTEP *timestep, double u[]);
void compute_dt_euler(GRID *grid, TIMESTEP *timestep, double D[]);
void solve(GRID *grid, double *A[], double b[]);
void iterative_solve(double X[], double *A[], double f[], int modes);
// pmain.c
void print_solution(char *string, GRID *grid, TIMESTEP *timestep, double u[]);
void initialize_states(GRID *grid, double u[]);
void choose_grid(GRID *grid);
void select_dt(TIMESTEP *timestep);

/*
*				pmain.c
*
*	Controls execution for the semiconductor process simulator.
*	The time integration method is TRBDF2.
*/

int main()
{
	double *u;		// current concentration array at time n
	GRID grid;
	TIMESTEP timestep;
	double t_old, t_max;
	int modes, step, MAX_STEPS;

	printf("\n\n1D semiconductor process simulation\n\n");

	timestep.t = 0.;
	fprintf(stderr, "Enter the max value of t in 1000 sec: ");
	scanf("%lf", &t_max);
	timestep.t_max = t_max;
	fprintf(stderr, "Enter the max number of timesteps: ");
	scanf("%d", &MAX_STEPS);
	fprintf(stderr, "Enter initial FACTOR & MAX_FACTOR for ");
	fprintf(stderr, "dt = FACTOR*dt_euler(t = 0): ");
	scanf("%lf %lf", &timestep.FACTOR, &timestep.MAX_FACTOR);
	fprintf(stderr, "Enter MIN_FACTOR for dt (e.g. 0.01): ");
	scanf("%lf", &timestep.MIN_FACTOR); // If dt<dt_min, program exits!

	printf("dt initial FACTOR = %g, MIN_FACTOR = %g, MAX_FACTOR = %g\n",
		timestep.FACTOR, timestep.MIN_FACTOR, timestep.MAX_FACTOR);

	choose_grid(&grid);
	modes = grid.modes = grid.Nx + 1;

	// allocate storage for array u
	u = (double *) alloc_vector(modes, DOUBLE);

	initialize_states(&grid, u);
	print_solution("initial data", &grid, &timestep, u);
	check_conservation(&grid, u);

			// timestepping loop
	timestep.step = 0;
	timestep.LAST_STEP = NO;
	timestep.REDO_STEP = NO;
	compute_first_dt(&grid, &timestep, u);
	printf("calculated timestep = %g\n", timestep.dt);
	for (step = 1; step <= MAX_STEPS && !timestep.LAST_STEP; step++) {
		timestep.step = step;
		do {
			t_old = timestep.t;
			select_dt(&timestep);
			printf("\ntimestep: %d, t = %g, ", step, timestep.t);
			printf("dt = %g, FACTOR = %g\n", timestep.dt,
				timestep.FACTOR);
			TRBDF2(&grid, &timestep, u);
			if (timestep.REDO_STEP) {
				timestep.t = t_old;
				if (timestep.FACTOR < timestep.MIN_FACTOR/2.) {
					printf("ERROR: dt < dt_min\n");
					exit(ERROR);
				}
			}
		} while (timestep.REDO_STEP);
		check_conservation(&grid, u);
	}
	print_solution("solution", &grid, &timestep, u);

	return 0;
}

void print_solution(char *string, GRID *grid, TIMESTEP *timestep, double u[])
{
	int Nx = grid->Nx, i;
	double x;
	static int first = YES;
	static FILE *file1, *file2;

	if (first) {
		first = NO;
		file1 = fopen("TRBDF2.sol0", "w");
		printf("\n%s at time t = %g in TRBDF2.sol0\n", string,
			timestep->t);
		for (i = 0; i <= Nx; i++) {
			x = grid->xmin + i*grid->dx;
			fprintf(file1, "%-12g %-12g\n", x, u[i]);
		}
		fclose(file1);
	}

	else {
		file2 = fopen("TRBDF2.sol", "w");
		printf("\n%s at time t = %g in TRBDF2.sol\n", string,
			timestep->t);
		for (i = 0; i <= Nx; i++) {
			x = grid->xmin + i*grid->dx;
			fprintf(file2, "%-12g %-12g\n", x, u[i]);
		}
		fclose(file2);
	}
}

void initialize_states(GRID *grid, double u[])
{
	int Nx = grid->Nx, i;
	double x, x0;

	// initialize problem
	x0 = grid->xmin + (grid->xmax-grid->xmin)/10.;
	for (i = 0; i <= Nx; i++) {
		x = grid->xmin + i*grid->dx;
		u[i] = 1.5*exp(-2.*sq(x-x0));
	}
}

void choose_grid(GRID *grid)
{
	double xmin, xmax;
	int Nx;

	fprintf(stderr, "Enter the min & max values of x in 0.1 microns: ");
	scanf("%lf %lf", &xmin, &xmax);
	grid->xmin = xmin;
	grid->xmax = xmax;
	fprintf(stderr, "Enter the number of dx: ");
	scanf("%d", &Nx);
	grid->Nx = Nx;
	printf("number of dx = %d\n", Nx);
	grid->dx = (xmax-xmin)/Nx;
}

void select_dt(TIMESTEP *timestep)
{
	double t_max = timestep->t_max;

	timestep->LAST_STEP = NO;
	if (timestep->t+timestep->dt > t_max) {
		timestep->dt = t_max - timestep->t;
		timestep->t = t_max;
		timestep->FACTOR = timestep->dt/timestep->dt_euler;
		timestep->LAST_STEP = YES;
	}
	else timestep->t += timestep->dt;
}

/*
*				parab.c
*
*	Contains the TRBDF2() parabolic solver & compute_dt().
*/

#define EPSILON		1.e-6			// for Newton's method
#define EPSILON2	0.5e-6			// for Newton's method
#define GAMMA		(2.-sqrt(2.))		// for TRBDF2

void TRBDF2(GRID *grid, TIMESTEP *timestep, double u[])
{
	int i, Nx = grid->Nx, modes = grid->modes;
	double dt;
	static double *u_mid;		// concentration array at n + gamma
	static double *u_new;		// concentration array at n + 1
	static double *F;		// current F array at time n
	static double *F_mid;		// F at n + gamma
	static double *F_new;		// F at n + 1
	static double *D;		// diffusion coefficient array
	static double *deltaD;		// (delta D)/(delta u)  array
	static double **J;		// Jacobian for TR & BDF2 steps
	static double *b;		// RHS for TR & BDF2 steps
	static int first = YES;

	if (first) {
		first = NO;
		u_mid = (double *) alloc_vector(modes, DOUBLE);
		u_new = (double *) alloc_vector(modes, DOUBLE);
		F = (double *) alloc_vector(modes, DOUBLE);
		F_mid = (double *) alloc_vector(modes, DOUBLE);
		F_new = (double *) alloc_vector(modes, DOUBLE);
		D = (double *) alloc_vector(modes, DOUBLE);
		deltaD = (double *) alloc_vector(modes, DOUBLE);
		J = (double **) alloc_matrix(modes, modes, DOUBLE);
		b = (double *) alloc_vector(modes, DOUBLE);
	}

	dt = timestep->dt;
	TR_step(GAMMA*dt, grid, timestep, u, u_mid, F, F_mid, D, deltaD, J, b);
	BDF2_step(dt, grid, timestep, u, u_mid, u_new, F_new, D, deltaD, J, b);
	compute_dt(grid, timestep, F, F_mid, F_new);
	if (!timestep->REDO_STEP) for (i = 0; i <= Nx; i++) u[i] = u_new[i];
}

// Note that the timestep is named dt_TR in the subroutine; called with
// argument GAMMA*dt above for TRBDF2
void TR_step(double dt_TR, GRID *grid, TIMESTEP *timestep, double u[],
	double u_mid[], double F[], double F_mid[], double D[],
	double deltaD[], double *J[], double b[])
{
	int i, step, Nx = grid->Nx, modes = grid->modes, MAX_NEWTONS = 8;
	double norm_change = 1., norm_residual = 1.;

	for (i = 0; i <= Nx; i++) u_mid[i] = u[i];
	compute_D(grid, u, D);
	compute_F(grid, u, D, F);

	printf("\nTR step:\n");
			// Newton iteration loop
	for (step = 1; step <= MAX_NEWTONS && norm_residual > EPSILON &&
			norm_change > EPSILON2; step++) {
		printf("\tNewton iteration: %d\n", step);

		compute_D(grid, u_mid, D);
		compute_F(grid, u_mid, D, F_mid);
		compute_deltaD(grid, u_mid, D, deltaD);
		compute_deltaF(grid, u_mid, D, deltaD, J);

		// evaluate jacobian
		compute_TR_J(dt_TR, grid, J);

		// evaluate residual
		compute_TR_b(dt_TR, grid, u, u_mid, F, F_mid, b);
		compute_norm(b, modes, &norm_residual);
		printf("\tnorm_residual in newton() = %g\n", norm_residual);

		solve(grid, J, b);

		// check for convergence
		compute_norm(b, modes, &norm_change);
		printf("\tnorm_change in newton() = %g\n", norm_change);

		update_states(grid, timestep, u_mid, b);
	}

	if (step > MAX_NEWTONS) {
		printf("\nERROR: Newton's method failed to converge after ");
		printf("%d iterations\n", MAX_NEWTONS);
		exit(ERROR);
	}
	if (norm_residual > EPSILON) {
		printf("\nWARNING: Newton's method converged on norm_change");
		printf(" = %g < %g,\nbut ", norm_change, EPSILON2);
		printf("norm_residual = %g > %g\n", norm_residual, EPSILON);
	}
}

void BDF2_step(double dt, GRID *grid, TIMESTEP *timestep, double u[],
	double u_mid[], double u_new[], double F_new[], double D[],
	double deltaD[], double *J[], double b[])
{
	int i, step, Nx = grid->Nx, modes = grid->modes, MAX_NEWTONS = 8;
	double norm_change = 1., norm_residual = 1.;

	for (i = 0; i <= Nx; i++) u_new[i] = u_mid[i];

	printf("\nBDF2 step:\n");
			// Newton iteration loop
	for (step = 1; step <= MAX_NEWTONS && norm_residual > EPSILON &&
			norm_change > EPSILON2; step++) {
		printf("\tNewton iteration: %d\n", step);

		compute_D(grid, u_new, D);
		compute_F(grid, u_new, D, F_new);
		compute_deltaD(grid, u_new, D, deltaD);
		compute_deltaF(grid, u_new, D, deltaD, J);

		// evaluate jacobian
		compute_BDF2_J(dt, grid, J);

		// evaluate residual
		compute_BDF2_b(dt, grid, u, u_mid, u_new, F_new, b);
		compute_norm(b, modes, &norm_residual);
		printf("\tnorm_residual in newton() = %g\n", norm_residual);

		solve(grid, J, b);

		// check for convergence
		compute_norm(b, modes, &norm_change);
		printf("\tnorm_change in newton() = %g\n", norm_change);

		update_states(grid, timestep, u_new, b);
	}

	if (step > MAX_NEWTONS) {
		printf("\nERROR: Newton's method failed to converge after ");
		printf("%d iterations\n", MAX_NEWTONS);
		exit(ERROR);
	}
	if (norm_residual > EPSILON) {
		printf("\nWARNING: Newton's method converged on norm_change");
		printf(" = %g < %g,\nbut ", norm_change, EPSILON2);
		printf("norm_residual = %g > %g\n", norm_residual, EPSILON);
	}
}

void compute_F(GRID *grid, double u[], double D[], double F[])
{
	int i, Nx = grid->Nx;
	double dx = grid->dx;

	for (i = 1; i < Nx; i++)
		F[i] = 0.5*( (D[i]+D[i+1])*(u[i+1]-u[i]) -
			(D[i-1]+D[i])*(u[i]-u[i-1]) )/sq(dx);
	// Neumann BCs
	F[0] = (D[0]+D[1])*(u[1]-u[0])/sq(dx);
	F[Nx] = (D[Nx]+D[Nx-1])*(u[Nx-1]-u[Nx])/sq(dx);
}

// delta_F*sq(dx)
void compute_deltaF(GRID *grid, double u[], double D[], double deltaD[],
	double *J[])
{
	int Nx = grid->Nx, i, j;

	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Nx; j++) J[i][j] = 0.;

	// Neumann BCs
	J[0][0] = -(D[1]+D[0])+deltaD[0]*(u[1]-u[0]);
	J[0][1] = (D[1]+D[0])+deltaD[1]*(u[1]-u[0]);
	J[Nx][Nx] = -(D[Nx-1]+D[Nx])+deltaD[Nx]*(u[Nx-1]-u[Nx]);
	J[Nx][Nx-1] = (D[Nx]+D[Nx-1])-deltaD[Nx-1]*(u[Nx]-u[Nx-1]);

	for (i = 1; i < Nx; i++) {
		J[i][i] = -0.5*(D[i+1]+2.*D[i]+D[i-1])+
			0.5*deltaD[i]*(u[i+1]-2.*u[i]+u[i-1]);
		J[i][i+1] = 0.5*(D[i+1]+D[i])+0.5*deltaD[i+1]*(u[i+1]-u[i]);
		J[i][i-1] = 0.5*(D[i]+D[i-1])-0.5*deltaD[i-1]*(u[i]-u[i-1]);
	}
}

void compute_TR_J(double dt_TR, GRID *grid, double *J[])
{
	int Nx = grid->Nx, i;
	double constant;

	// calculate (1 - Jdt/2)
	constant = 0.5*dt_TR/sq(grid->dx);
	for (i = 0; i <= Nx; i++) {
		J[i][i] = 1. - constant*J[i][i];
		if (i > 0) J[i][i-1] *= -constant;
		if (i < Nx) J[i][i+1] *= -constant;
	}
}

void compute_BDF2_J(double dt, GRID *grid, double *J[])
{
	int Nx = grid->Nx, i;
	double constant;

	// calculate (1 - dt*(1-GAMMA)*J/(2-GAMMA))
	constant = dt*(1.-GAMMA)/((2.-GAMMA)*sq(grid->dx));
	for (i = 0; i <= Nx; i++) {
		J[i][i] = 1. - constant*J[i][i];
		if (i > 0) J[i][i-1] *= -constant;
		if (i < Nx) J[i][i+1] *= -constant;
	}
}

void compute_TR_b(double dt_TR, GRID *grid, double u[], double u_mid[],
	double F[], double F_mid[], double b[])
{
	int Nx = grid->Nx, i;

	for (i = 0; i <= Nx; i++)
		b[i] = -(u_mid[i]-u[i]) + 0.5*dt_TR*(F_mid[i]+F[i]);
}

void compute_BDF2_b(double dt, GRID *grid, double u[], double u_mid[],
	double u_new[], double F_new[], double b[])
{
	int Nx = grid->Nx, i;

	for (i = 0; i <= Nx; i++)
		b[i] = -(u_new[i]-u_mid[i]/(GAMMA*(2.-GAMMA))+
			sq(1.-GAMMA)*u[i]/(GAMMA*(2.-GAMMA))) +
			dt*(1.-GAMMA)*F_new[i]/(2.-GAMMA);
}

void update_states(GRID *grid, TIMESTEP *timestep, double u[], double b[])
{
	int modes = grid->modes, i;

	for (i = 0; i < modes; i++) u[i] += b[i];
	compute_norm(u, modes, &timestep->norm_u);
}

#define DEBUG_dt TRUE
void compute_dt(GRID *grid, TIMESTEP *timestep, double F[], double F_mid[],
	double F_new[])
{
	int modes = grid->modes, i;
	double e_l; // = ||e_l||
	double dt = timestep->dt, dt_new;
	double dt_max = timestep->dt_max;
	double r, k = (-3.*sq(GAMMA)+4.*GAMMA-2.)/(12.*(2.-GAMMA));
	double EPSILON_REL = 1.e-6, EPSILON_ABS = 1.e-12;

	e_l = 0.;
	for (i = 0; i < modes; i++)
		e_l += fabs(2.*k*dt*(F[i]/GAMMA - F_mid[i]/(GAMMA*(1.-GAMMA)) +
			F_new[i]/(1.-GAMMA)));
	e_l /= modes;
	r = e_l/(EPSILON_REL*timestep->norm_u + EPSILON_ABS);
	if (DEBUG_dt) printf("TRBDF2 r = %g\n", r);

	if (r <= 2.) {
		timestep->REDO_STEP = NO;
		dt_new = dt/pow(r, 1./3.);
		dt = min(dt_new, 5.*dt);
	}
	else {
		timestep->REDO_STEP = YES;
		dt /= 2.;
		printf("WARNING: TRBDF2 r = %g > 2\n", r);
		printf("dt /= 2.; redoing timestep %d\n", timestep->step);
	}
	dt = min(dt, dt_max);

	timestep->dt = dt;
	timestep->FACTOR = dt/timestep->dt_euler;
}

void compute_first_dt(GRID *grid, TIMESTEP *timestep, double u[])
{
	int modes = grid->modes;
	double *D;

	D = (double *) alloc_vector(modes, DOUBLE);
	compute_D(grid, u, D);
	compute_dt_euler(grid, timestep, D);
	free((char *) D);
	timestep->dt = timestep->FACTOR*timestep->dt_euler;
	timestep->dt_min = timestep->MIN_FACTOR*timestep->dt_euler;
	timestep->dt_max = timestep->MAX_FACTOR*timestep->dt_euler;
}

// maximum dt for forward Euler
void compute_dt_euler(GRID *grid, TIMESTEP *timestep, double D[])
{
	int i, Nx = grid->Nx;
	double dt_euler, max_D;
	double dx = grid->dx;

	max_D = 0.;

	for (i = 0; i <= Nx; i++) max_D = max(D[i], max_D);
	dt_euler = sq(dx)/(2.*max_D);

	timestep->dt_euler = dt_euler;
}

void solve(GRID *grid, double *A[], double b[])
{
	int i, modes = grid->modes;
	static double *X;
	static int first = YES;

	// ITERATIVE SOLVE
	if (first) {
		first = NO;
		X = (double *) alloc_vector(modes, DOUBLE);
	}
	iterative_solve(X, A, b, modes);
	for (i = 0; i < modes; i++) b[i] = X[i];
}

void iterative_solve(double X[], double *A[], double f[], int modes)
// A X = f
#define N_iterations	1000 // maximum number of iterations
{
	double residual = 1., residual0, sum;
	double EPSILON_GS = 1.e-6;
	int i, j, step;

			// Gauss-Seidel iteration
	residual0 = 0.;
	for (i = 0; i < modes; i++) residual0 += fabs(f[i]);
	residual0 /= modes;

	for (step = 1; step <= N_iterations && residual > EPSILON_GS*residual0;
			step++) {
		for (i = 0; i < modes; i++) {
			X[i] = f[i];
			for (j = 0; j < modes; j++)
				if (j != i) X[i] += -A[i][j]*X[j];
			X[i] /= A[i][i];
		}

		residual = 0.;
		for (i = 0; i < modes; i++) {
			sum = -f[i];
			for (j = 0; j < modes; j++) sum += A[i][j]*X[j];
			residual += fabs(sum);
		}
		residual /= modes;
	}
	printf("Gauss-Seidel iterations = %d\n", --step);
	if (residual > EPSILON_GS*residual0)
		printf("ERROR: Gauss-Seidel residual r/r0= %g > EPSILON_GS\n",
			residual);
}

/*
*				diffuse.c
*
*	Contains the phenomenological diffusion model for BORON for
*	process simulations.
*
*	The diffusion coefficient is modelled as
*		D_phys = alpha_phys*u_phys + beta_phys (in physical units).
*	In computational units,
*		D = D_phys/D_scale
*		  = (alpha*u_scale)*u + beta
*	For example, if
*		alpha_phys = 1.5*10^-33 cm^5/sec &
*		beta_phys = 2*10^-15 cm^2/sec, then
*	then
*		alpha = alpha_phys*u_scale/D_scale
*		      = 1.5
*		beta = beta_phys/D_scale
*		     = 2*10^-2
*	where u_scale & D_scale are defined in pdecs.h as
*	#define u_scale		1.e20
*	#define D_scale		1.e-13
*/

double alpha = 1.5, beta = 0.02;

void compute_D(GRID *grid, double u[], double D[])
{
	int i, modes = grid->modes;

	for (i = 0; i < modes; i++) D[i] = alpha*u[i] + beta;
}

void compute_deltaD(GRID *grid, double u[], double D[], double deltaD[])
{
	int i, modes = grid->modes;

	for (i = 0; i < modes; i++) deltaD[i] = alpha;
}

/*
*				cons.c
*
*	Checks the conservation of impurities.
*/

void check_conservation(GRID *grid, double u[])
// For 100 dx & 100 dt, (original-final)/original = 2 parts in 10^11
{
	double total;
	static double original_total;
	static int first = YES;

	total = compute_total(grid, u);

	if (first) {
		first = NO;
		printf("total u = %g\n", original_total = total);
	}
	else printf("total u = %g\t(original-current)/original = %g\n",
		total, (original_total-total)/original_total);
}

double compute_total(GRID *grid, double u[])
{
	double dx = grid->dx, sum;
	int Nx = grid->Nx, i;

	sum = 0.;
	for (i = 0; i < Nx; i++) sum += u[i] + u[i+1];
	sum *= 0.5*dx;

	return sum;
}

/*
*				Library Functions
*	Copyright 2015 by Carl L. Gardner.  All rights reserved.
*/

void compute_norm(double v[], int N, double *norm)
{
	int i;

	*norm = 0.;
	for (i = 0; i < N; i++) *norm += fabs(v[i]);
	*norm /= N;
}

POINTER Alloc(unsigned N_bytes)
{
	POINTER p;

	p = malloc(N_bytes);

	if (p == (POINTER) NULL) {
		perror("ERROR in Alloc(): malloc() returned NULL pointer");
		exit(ERROR);
	}
	return p;
}

POINTER alloc_vector(int N, unsigned element_size)
{
	return Alloc((unsigned) (N*element_size));
}

POINTER alloc_matrix(int N_rows, int N_columns,	unsigned element_size)
{
	int i, space, N_pointers;
	POINTER pA, array_origin;

		// due to alignment requirement
	N_pointers = (N_rows%2) ? N_rows+1 : N_rows;

	space = N_pointers*sizeof(POINTER) + N_rows*N_columns*element_size;
	pA = Alloc((unsigned) space);

	array_origin = ((char *) (pA)) + N_pointers*sizeof(POINTER);

	for (i = 0; i < N_rows; i++)
		((POINTER *) pA)[i] = ((char *) array_origin) +
			i*N_columns*element_size;

	return pA;
}
