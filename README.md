# statistical-estimation-theory

Archival MATLAB code from ASE 381P: Statistical Estimation Theory at the University of Texas at Austin with Dr. Todd Humphreys.

I originally took this course in Fall 2019 and reorganized/published my code in early 2022. The course had a lot of derivation (not included here), but problems with MATLAB computation are included here. Most of this is problem-specific, but should give starter code for a range of estimation techniques including:
- Detection Theory/Hypothesis testing
- Batch Least-Squares Estimation
- Kalman Filtering
- Unscented Kalman Filtering
- Information Filtering
- Particle Filtering

Note: Course textbook was Bar-Shalom *et al*, "Estimation with Applications to Tracking and Navigation: Theory, Algorithms, and Software"

# Organization
This repo is organized as follows:
- `/problem-sets` and `exams` contain the problem sets & exams from the course in .PDF format.
- `/data` contains data for use in the problems, such as simulated measurements.
- `/functions` contains reusable helper functions.
- `/scripts` contains the caller scripts for the problems in the problem sets and exams. A brief description of each MATLAB problem is given in the following section.

# Problem descriptions (`/scripts`)
- **ps1_p8.m**: Central limit theorem/random sampling problem. As you average uniformly sampled numbers, the resulting distribution approximates a Gaussian.
- **ps2_p3.m**: Matrix conditioning exercise. Compares differences in numerical precision (matrix condition) between analytical, left division, and least squares solutions of linear systems.
- **ps2_p5.m**: Hypothesis testing using Chi-squared distributions. Finds probability of differentiating one zero-mean Gaussian-distributed variable from another.
- **ps2_bs_21.m**: Problem 2-1 from Bar-Shalom which compares minimum mean-square error (MMSE) and maximum-a-priori (MAP) results for a simple 1D Bayesian estimator. 
- **ps3_p2.m**: Comparison of results for a batch Least Squares estimator with a Square Root Information Based Least Squares (SRIBLS) method. z = measurement vector, R = measurement noise covariance, and H = linearized measurement sensitivity matrix (Z=Hx). SRIBLS avoids a complicated matrix inversion by performing a Cholesky decomposition on the measurement noise (R) matrix, and instead performs the inversion at the end to recover the state estimate x_hat. Gives similar results but square root works better when matrices are poorly conditioned (i.e. near singular). 


# Known errors/needed fixes
- ps2_bs_21.m needs reworked. Gives incorrect numbers and could probably use a plot or two to really illustrate the point of the problem. Parts of it could potentially be turned into a function for reusability.
- 