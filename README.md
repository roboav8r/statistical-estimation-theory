# statistical-estimation-theory

Archival MATLAB code from ASE 381P: Statistical Estimation Theory at the University of Texas at Austin with Dr. Todd Humphreys.

I originally took this course in Fall 2019 and reorganized/published my code in early 2022. The course had a lot of derivation (not included here), but any problems with MATLAB computation are included. Most of this is problem-specific and designed to illustrate a point about probability/estimation techniques instead of usable code (most programming languages have libraries for this nowadays...). However, it can still be used as starter code/idea generation for a range of estimation techniques including:
- Detection Theory/Hypothesis testing
- Nonlinear equation solutions (Newton's method, Gauss-Newton)
- Batch Least-Squares Estimation
- Kalman Filtering
- Unscented Kalman Filtering
- Information Filtering
- Smoothing
- Particle Filtering

Note: Course textbook was Bar-Shalom *et al*, "Estimation with Applications to Tracking and Navigation: Theory, Algorithms, and Software"

# Organization
This repo is organized as follows:
- `/problem-sets` and `/exams` contain the problem sets & exams from the course in .PDF format.
- `/data` contains data for use in the problems, such as simulated measurements.
- `/functions` contains reusable helper functions.
- `/scripts` contains my caller scripts for the problems in the problem sets and exams. A brief description of each MATLAB problem is given in the following section.

# Problem descriptions (`/scripts`)
MATLAB scripts and a basic description of their problem are listed below. The naming convention is that ps5_p2.m refers to problem set 5, problem 2.

## Problem set 1: Probability
- **ps1_p8.m**: Central limit theorem/random sampling problem. As you average a higher number of uniformly sampled numbers, the resulting distribution approximates a Gaussian.

## Problem set 2: Hypothesis testing and Bayesian estimators
- **ps2_p3.m**: Matrix conditioning exercise. Compares differences in numerical precision (matrix condition) between analytical, left division, and least squares solutions of linear systems.
- **ps2_p5.m**: Hypothesis testing using Chi-squared distributions. Finds probability of differentiating one zero-mean Gaussian-distributed variable from another.
- **ps2_bs_21.m**: Problem 2-1 from Bar-Shalom which compares minimum mean-square error (MMSE) and maximum-a-priori (MAP) results for a simple 1D Bayesian estimator. 

## Problem sets 3 & 4: Nonlinear Least Squares (NLLS) estimators / static parameter estimation 
- **ps3_p2.m**: Comparison of results for a set of linear equations solved with two method: overconstrained matrix algebra, and with a Square Root Information Based Least Squares (SRIBLS) method. $z$ = measurement vector, $R$ = measurement noise covariance, and $H$ = linearized measurement sensitivity matrix ($Z=Hx$). SRIBLS avoids a matrix inversion by performing a Cholesky decomposition on the measurement noise matrix $R$, and instead performs the inversion at the end to recover the state estimate $\hat{x}$. Gives similar results but square root works better when matrices are poorly conditioned (i.e. near singular). 
- **ps4_p1.m**: Newton's method for solving a set of nonlinear equations. Different initial guesses result in different solutions.
- **ps4_p3.m**: Gauss-Newton method for finding the parameters (amplitude, frequency, phase) of a sinusoidal signal. Uses a least-squares cost function.
- **ps4_p4.m**: Estimation of missile ballistic trajectory parameters from two radar ground station range measurements. Uses Gauss-Newton and least-squares cost function.
- **ps4_p5.m**: Visualization of a Markov process function. Helper function for a derivation problem.
- **ps4_p6.m**: Effects of step-size and in NLLS cost reduction. Shows how step size can result in convergence or divergence of error depending on how it's chosen.

## Problem sets 5 & 6: Kalman Filtering, Square Root Information Filtering, and Smoothing
- **ps5_p3.m**: Kalman Filtering of a linear time invariant (LTI) system with known state, measurement, and noise matrices. Plots estimated state variable values & standard deviation/variance over time.
- **ps5_p4.m**: Same system as problem 3, but compares my computed values to MATLAB's `kalman.m`-computed values. Also checks system stability/convergence by analyzing the Eigenvalues of the error transition matrix and ensuring they have magnitude less than 1.
- **ps5_p5.m**: Evaluating three different LTI Kalman Filters with different process noise $Q$ and measurement noise $R$ values where the true parameters are unknown. Uses Chi-squared distributions and computes the average $\epsilon_\nu$ to determine the correct model.
- **ps5_p6.m**: Evaluating the consistency of a LTI KF. Filters are evaluated on synthetic measurement data generated via Monte Carlo simulation.
- **ps6_p1.m**: Evaluating the consistency (like ps5_p6) of a filter that has known inconsistent parameters.
- **ps6_p2.m**: Comparing Kalman Filter with Square Root Information Filter (SRIF) for the same linear system. Like `ps3_p2`, the SRIF performs better when matrices are not well-conditioned.
- **ps6_p3.m**: Smoothing problem applied to the same system in `ps6_p2`. Uses Kalman Filtering, Kalman Smoothing, Information Filtering, and Information Smoothing. 

## Problem set 7: Nonlinear Estimation & the Extended Kalman Filter
- **ps7_p1.m**: Comparison of analyical and linearized orbit propagation models. The linearized system uses MATLAB's `ode45` propagator. The difference  in the analytical and propagated orbit radii is less than 1 percent.

## Problem set 8: Multiple-Model Kalman Filters, Sigma-Point (Unscented) Filter, and Particle Filtering
- **ps8_prob_II.m**: Sigma-Point/Unscented Kalman filter problem with comparison to EKF.

# Known errors/needed fixes

Here are a few things I noticed on the 2022 review that I'd like to correct for future usage:

- Overall generalization and modularity of common functions.
- Make constants and parameters global variables instead of function inputs.
- Use consistent variable names across the board (e.g. sometimes Kalman gain is $W_ss$, sometimes it's just $W$ or $K$)
- ps2_p5.m (Simple hypothesis test for two zero-mean gaussian variables with different standard deviations): Make a function that takes $\sigma_0$, $\sigma_1$, and probability of detection $P_D$ as inputs, then returns the $\sigma_D$ decision value that gives the desired detection rate, and have this be the simple heuristic ($\lt \sigma_D$ is classified as $\sigma_0$'s variable, $\gt \sigma_D$ is classified as $\sigma_1$'s variable). Also plots would work well here to show the point of the problem: vertical line at decision/hypothesis test sigma value.
- ps2_bs_21.m (MAP vs. MMSE estimator) needs reworked. Gives incorrect numbers and could probably use a plot or two to really illustrate the point of the problem. Parts of it could potentially be turned into a function for reusability.
- ps4_p4.m: Add initial measurements to charts. Label charts. 
- Make generalized Gauss-Newton solver and matrix normalization helper functions. Use them for ps4_p3 and ps4_p4.
- ps4_p6: Implement variable step size (alpha).
- ps5_p5 and ps5_p6: Revisit these and clearly define the problem and objectives. I don't think the Monte Carlo averaging was implemented as intended.
- `ps8_p1.m`, a multiple Kalman Filter problem, and `ps8_prob_III.m`, a particle filtering problem, are missing. Find and add to repo.
