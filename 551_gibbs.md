STATS 551 - Posterior Computation Techniques: Discretization & Gibbs
Sampler
================
GSI: Bach Viet Do and Vincenzo LoffredoProfessor: Long Nguyen
06 March, 2021

# Motivation

Le ![\\theta](https://latex.codecogs.com/png.latex?%5Ctheta "\theta") be
a continous random variable with the prior
![p(\\theta)](https://latex.codecogs.com/png.latex?p%28%5Ctheta%29 "p(\theta)")
and likelihood
![p(x \| \\theta)](https://latex.codecogs.com/png.latex?p%28x%20%7C%20%5Ctheta%29 "p(x | \theta)"),
the goal of Bayesian inference is to compute the posterior:

-   If we have conjugacy, for example,
    ![p(x \| \\theta) \\sim \\text{Bernoulli}(x \| \\theta)](https://latex.codecogs.com/png.latex?p%28x%20%7C%20%5Ctheta%29%20%5Csim%20%5Ctext%7BBernoulli%7D%28x%20%7C%20%5Ctheta%29 "p(x | \theta) \sim \text{Bernoulli}(x | \theta)")
    and
    ![p(\\theta) \\sim \\text{Beta}(\\theta \| a, b)](https://latex.codecogs.com/png.latex?p%28%5Ctheta%29%20%5Csim%20%5Ctext%7BBeta%7D%28%5Ctheta%20%7C%20a%2C%20b%29 "p(\theta) \sim \text{Beta}(\theta | a, b)")
    then
    ![p(\\theta \| x) = \\cfrac{\\theta^{x + a - 1}(1-\\theta)^{b + 1 - x -1}}{\\int p(x,\\theta) d\\theta}](https://latex.codecogs.com/png.latex?p%28%5Ctheta%20%7C%20x%29%20%3D%20%5Ccfrac%7B%5Ctheta%5E%7Bx%20%2B%20a%20-%201%7D%281-%5Ctheta%29%5E%7Bb%20%2B%201%20-%20x%20-1%7D%7D%7B%5Cint%20p%28x%2C%5Ctheta%29%20d%5Ctheta%7D "p(\theta | x) = \cfrac{\theta^{x + a - 1}(1-\theta)^{b + 1 - x -1}}{\int p(x,\theta) d\theta}").
    We recognize the signature of Beta distribution in the numerator and
    because the Bayes rule ensures a proper distribution, we conclude
    that
    ![\\int p(x, \\theta) dx = \\cfrac{\\Gamma(x + a)\\Gamma(b + 1 - x)}{\\Gamma(a + b \_ 1)}](https://latex.codecogs.com/png.latex?%5Cint%20p%28x%2C%20%5Ctheta%29%20dx%20%3D%20%5Ccfrac%7B%5CGamma%28x%20%2B%20a%29%5CGamma%28b%20%2B%201%20-%20x%29%7D%7B%5CGamma%28a%20%2B%20b%20_%201%29%7D "\int p(x, \theta) dx = \cfrac{\Gamma(x + a)\Gamma(b + 1 - x)}{\Gamma(a + b _ 1)}"),
    and so bypass computing the integration.

-   On the other hand, let’s say
    ![p(x \| \\theta) = \\cfrac{1}{2}\\cdot \\text{N}(x \| \\theta, 1) + \\cfrac{1}{2}\\cdot \\text{N}(x \| -\\theta, 1)](https://latex.codecogs.com/png.latex?p%28x%20%7C%20%5Ctheta%29%20%3D%20%5Ccfrac%7B1%7D%7B2%7D%5Ccdot%20%5Ctext%7BN%7D%28x%20%7C%20%5Ctheta%2C%201%29%20%2B%20%5Ccfrac%7B1%7D%7B2%7D%5Ccdot%20%5Ctext%7BN%7D%28x%20%7C%20-%5Ctheta%2C%201%29 "p(x | \theta) = \cfrac{1}{2}\cdot \text{N}(x | \theta, 1) + \cfrac{1}{2}\cdot \text{N}(x | -\theta, 1)")
    and
    ![p(\\theta) \\sim N(0,1)](https://latex.codecogs.com/png.latex?p%28%5Ctheta%29%20%5Csim%20N%280%2C1%29 "p(\theta) \sim N(0,1)"),
    then
    ![p(\\theta \| x) = \\cfrac{\\cfrac{1}{4\\pi} \\cdot \\exp\\left( - \\cfrac{(x - \\theta)^2 + \\theta^2}{2}\\right) + \\cfrac{1}{4\\pi} \\cdot \\exp\\left( - \\cfrac{(x + \\theta)^2 + \\theta^2}{2}\\right)}{\\int p(x, \\theta) d\\theta}](https://latex.codecogs.com/png.latex?p%28%5Ctheta%20%7C%20x%29%20%3D%20%5Ccfrac%7B%5Ccfrac%7B1%7D%7B4%5Cpi%7D%20%5Ccdot%20%5Cexp%5Cleft%28%20-%20%5Ccfrac%7B%28x%20-%20%5Ctheta%29%5E2%20%2B%20%5Ctheta%5E2%7D%7B2%7D%5Cright%29%20%2B%20%5Ccfrac%7B1%7D%7B4%5Cpi%7D%20%5Ccdot%20%5Cexp%5Cleft%28%20-%20%5Ccfrac%7B%28x%20%2B%20%5Ctheta%29%5E2%20%2B%20%5Ctheta%5E2%7D%7B2%7D%5Cright%29%7D%7B%5Cint%20p%28x%2C%20%5Ctheta%29%20d%5Ctheta%7D "p(\theta | x) = \cfrac{\cfrac{1}{4\pi} \cdot \exp\left( - \cfrac{(x - \theta)^2 + \theta^2}{2}\right) + \cfrac{1}{4\pi} \cdot \exp\left( - \cfrac{(x + \theta)^2 + \theta^2}{2}\right)}{\int p(x, \theta) d\theta}").
    It’s now not obvious what kind of distribution the numerator is and
    so the integration
    ![\\int p(x,\\theta) d\\theta](https://latex.codecogs.com/png.latex?%5Cint%20p%28x%2C%5Ctheta%29%20d%5Ctheta "\int p(x,\theta) d\theta")
    is no longer immediate.

As we can not work out the exact form the posterior, how do we sample or
perform computation with the posterior ? Today, we discuss two general
techniques to solve this problem: Posterior Discretization and Gibbs
Sampler/MCMC.

# Posterior Discretization

The discretization idea is to “pretend” that
![\\theta](https://latex.codecogs.com/png.latex?%5Ctheta "\theta") is
discrete and assumes a finite set of values (a “grid”)
![\\{\\theta\_1, \\theta\_2, \\ldots, \\theta\_S \\}](https://latex.codecogs.com/png.latex?%5C%7B%5Ctheta_1%2C%20%5Ctheta_2%2C%20%5Cldots%2C%20%5Ctheta_S%20%5C%7D "\{\theta_1, \theta_2, \ldots, \theta_S \}"),
and so the integration can be reduced into a finite sum. Our hope is by
choosing an approriate large S and reasonable range of the grid
![\\int p(x, \\theta) d\\theta \\approx \\sum\_{s=1}^S p(x, \\theta\_s) = \\sum\_{s=1}^S p(\\theta\_s) p(x \| \\theta\_s)](https://latex.codecogs.com/png.latex?%5Cint%20p%28x%2C%20%5Ctheta%29%20d%5Ctheta%20%5Capprox%20%5Csum_%7Bs%3D1%7D%5ES%20p%28x%2C%20%5Ctheta_s%29%20%3D%20%5Csum_%7Bs%3D1%7D%5ES%20p%28%5Ctheta_s%29%20p%28x%20%7C%20%5Ctheta_s%29 "\int p(x, \theta) d\theta \approx \sum_{s=1}^S p(x, \theta_s) = \sum_{s=1}^S p(\theta_s) p(x | \theta_s)")
and
![p(\\theta = \\theta' \| x) = \\cfrac{p(\\theta')p(x \| \\theta')}{\\sum\_{s=1}^S p(\\theta\_s) p(x \| \\theta\_s)}](https://latex.codecogs.com/png.latex?p%28%5Ctheta%20%3D%20%5Ctheta%27%20%7C%20x%29%20%3D%20%5Ccfrac%7Bp%28%5Ctheta%27%29p%28x%20%7C%20%5Ctheta%27%29%7D%7B%5Csum_%7Bs%3D1%7D%5ES%20p%28%5Ctheta_s%29%20p%28x%20%7C%20%5Ctheta_s%29%7D "p(\theta = \theta' | x) = \cfrac{p(\theta')p(x | \theta')}{\sum_{s=1}^S p(\theta_s) p(x | \theta_s)}").

## Example 1

In this example, we’ll sample from the posterior of a mixture of
Gaussian distribution. Assume that the data sample is generated as
follows:

![X\_1,\\ldots, X\_n \\sim 0.5 \\cdot \\text{N(-5, 1)} + 0.5 \\cdot \\text{N(5,1)} ](https://latex.codecogs.com/png.latex?X_1%2C%5Cldots%2C%20X_n%20%5Csim%200.5%20%5Ccdot%20%5Ctext%7BN%28-5%2C%201%29%7D%20%2B%200.5%20%5Ccdot%20%5Ctext%7BN%285%2C1%29%7D%20 "X_1,\ldots, X_n \sim 0.5 \cdot \text{N(-5, 1)} + 0.5 \cdot \text{N(5,1)} ")

We define the following model for the data:

Recall from Vincenzo’s lab, we generate this type of data by introducing
a (latent) variable z.

``` r
generate_mixture_model <- function(n = 1000, mu0 = 0, pi0 = 0.5) {
  x = rep(0, n) # initialize a sample of size n
  mu = c(-mu0, mu0)
  
  for(i in 1:n) {
    zi = rbinom(1, 1, pi0) + 1 # z[i] specifies which cluster x[i] belongs to
    x[i] = rnorm(1, mu[zi], 1)
  }
  
  x
}

x =  generate_mixture_model(n = 1000, mu0 = 5, pi0 = 0.5)

require(ggplot2)
ggplot(data.frame(x = x), aes(x=x)) + geom_histogram(color="black", fill="white", binwidth = 0.25) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Now, give a sample
![x\_1, \\ldots, x\_n](https://latex.codecogs.com/png.latex?x_1%2C%20%5Cldots%2C%20x_n "x_1, \ldots, x_n"),
we are interested in
![p(\\theta \| x\_1, \\ldots, x\_n)](https://latex.codecogs.com/png.latex?p%28%5Ctheta%20%7C%20x_1%2C%20%5Cldots%2C%20x_n%29 "p(\theta | x_1, \ldots, x_n)").
Let’s apply the disretization idea:

``` r
compute_post_theta_pdf <- function(theta, x, grid) {
  n = length(x)
  # calculate the log of  numerator
  log_up = dnorm(theta, mean = 10, log = T) +  sum(log(0.5 * dnorm(x, mean = theta) + 0.5 * dnorm(x, mean = -theta)))
  
  # approximating the integration with a finite sum of prior * likelihood
  # at each grid point
  down = 0
  for (i in 1:length(grid)) {
    log_ll =  sum(log(0.5 * dnorm(x, mean = grid[i]) + 0.5 * dnorm(x, mean = -grid[i]))) # add log likelihood at mean = theta_s
    log_joint = dnorm(grid[i], mean = 10, log = T) + log_ll # add prior for theta_s
    down = down + exp(log_joint - log_up)
  }
  
  1.00 / down
}
```

``` r
# compute posterior at a range of values
grid = seq(-10, 10, 0.1)
post_grid = rep(0, length(grid))

for(i in 1:length(grid)) {
  post_grid[i] = compute_post_theta_pdf(grid[i], x, grid)
}

# plot out the posterior density
ggplot(data = data.frame(x = grid, y = post_grid)) + geom_line(aes(x,y)) + xlab("theta") + ylab("p(theta | x)") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

**Remark** : From the posterior plot, the posterior density concentrates
around the value
![\\theta = 5](https://latex.codecogs.com/png.latex?%5Ctheta%20%3D%205 "\theta = 5")
which is the ground truth.

## Example 2

Recall from the semi-conjugate example in the note section 6.2 where we
derive a Gibbs Sampler to sample from the posterior
![p(\\theta, \\tau \| y)](https://latex.codecogs.com/png.latex?p%28%5Ctheta%2C%20%5Ctau%20%7C%20y%29 "p(\theta, \tau | y)").
In this lab, let’s redo it with the discrete approximation technique.
Denote
![\\tau = \\cfrac{1}{\\sigma^2}](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%20%5Ccfrac%7B1%7D%7B%5Csigma%5E2%7D "\tau = \cfrac{1}{\sigma^2}").
The model is defined as:

``` r
compute_joint_post_pdf <- function(y, mean.grid, prec.grid) {
  mu0 = 1.9 ; t20 = 0.95^2 ; s20 = .01 ; nu0 = 1
  
  post.grid = matrix(nrow=G, ncol = H)
  
  for(g in 1:G){
    for(h in 1:H) {
      post.grid[g, h] = dnorm (mean.grid[g] , mu0, sqrt(t20)) * dgamma(
        prec.grid[h],  nu0 / 2 , s20 * nu0 /2 ) * prod(dnorm( y, mean.grid[g] , 1/sqrt(prec.grid[h])))
    }
  }
  
  post.grid
}

y= c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)

ggplot(data.frame(x = y), aes(x=x)) + geom_histogram(color="black", fill="white", binwidth = 0.05) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
G = 100; H = 100
mean.grid = seq(1.505, 2.00, length = G)
prec.grid = seq(1.75, 175, length = H)

# compute the joint posterior (mu, tau) over a range of values
post.grid = compute_joint_post_pdf(y, mean.grid, prec.grid )
```

``` r
# plot the contour of the joint posterior (mu, tau) | y
require(reshape2)
post_contour = melt(post.grid)
colnames(post_contour) = c("theta", "tau", "post")
ggplot(post_contour, aes(x = theta, y = tau, z = post)) +
         stat_contour(geom = "polygon", aes(fill = ..level..)) +
         geom_tile(aes(fill = post)) +
         stat_contour(bins = 15) + ggtitle("Posterior (theta,tau) | y")  + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# plot the marginal posterior theta | y
post_theta = apply(post.grid, 1, sum)
ggplot(data = data.frame(x =mean.grid, y = post_theta)) + geom_line(aes(x,y)) + xlab("theta") + ylab("p(theta | y)") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# plot the marginal posteior tau | y
post_tau = apply(post.grid, 2, sum)
ggplot(data = data.frame(x = prec.grid, y = post_tau)) + geom_line(aes(x,y)) + xlab("tau") + ylab("p(tau | y)") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Gibbs Sampler

Another idea to sample from an intractable posterior distribution is to
create a Markov Chain of which the posterior is the stationary
distribution. Markov Chain theory guarantees that the chain converges to
the stationary distribution as the number of iterations goes to
infinity. Gibbs Sampler is a special case of Monte-Carlo Markov Chain
and covered in depth in the last two lectures. Today, we will go through
some concrete examples with coding.

From my experience, I recommend to code your Gibbs Sampler as an R OOP
object The advantage of this style of coding is that you can save or
checkpoint your Gibbs states to prevent losing everything in the event
of crashing. Moreover, you can also resume running later if there are
evidence of non-convergence. Lastly, it’s also easy to move the saved
Gibbs Sampler among diffeernt machines. Let’s first take a look into
Object Oriented Programming in R.

## Example 3

Let’s see how we use Gibbs Sampler for the problem in Example 2. Denote
![\\tau = \\cfrac{1}{\\sigma^2}](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%20%5Ccfrac%7B1%7D%7B%5Csigma%5E2%7D "\tau = \cfrac{1}{\sigma^2}").
Recall from the lecture note,

where

The implementation for GIbbs Sampler is as follows:

``` r
set.seed(1)

y= c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08) # data
y_bar = mean(y); y_var = var(y); n = length(y) # empirical mean and variance

mu0 = 1.9 ; t20 = 0.95^2 ; s20 = .01 ; nu0 = 1 # set hyperparmeter

sample_mu <- function() {
  taun = 1.0 / (1.0/t20 + n * tau)
  mun = taun * ( mu0 /t20 + n * y_bar * tau)
  
  mu <<- rnorm(1, mun, sqrt(taun))
}

sample_tau = function() {
  nun = nu0 + n
  s2n = (nu0 * s20 + (n - 1) * y_var  + n * (y_bar - mu)^2)/ nun
  
  tau <<- rgamma(1, nun/2, nun * s2n / 2)
}

gibbs_iter = function() {
  sample_mu()
  sample_tau()
}
```

``` r
niters = 1000

mu = y_bar; tau = 1 / y_var # set intial values

for(i in 1:niters) {
  gibbs_iter() # run Gibbs Sampler until convergence
}
```

``` r
N = 1000
x = matrix(nrow = N, ncol = 2) # create a sample of size N

# extract sample from conveged Gibbs
for(n in 1:N) {
  gibbs_iter()
  x[n, 1] = mu
  x[n, 2] = tau
}
```

``` r
require(ggplot2)
ggplot(data.frame(mu = x[,1]), aes(x=mu)) + geom_histogram(color="black", fill="white", binwidth = 0.005) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggplot(data.frame(tau = x[,2]), aes(x=tau)) + geom_histogram(color="black", fill="white", binwidth = 1) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

    ## [1] "95% CI for mu is :"

    ##     2.5%      50%    97.5% 
    ## 1.786751 1.804634 1.822816

    ## [1] "95% CI for tau is :"

    ##     2.5%      50%    97.5% 
    ## 47.11475 59.48439 71.22411

    ## [1] "95% CI for sigma.2 is :"

    ##       2.5%        50%      97.5% 
    ## 0.01404020 0.01681113 0.02122487

## Example 4

In the this example, let’s see how we can use Gibbs Sampler to sample
from a bivarite normal distribution. A bivariate normal distribution
with mean vector
![\\mu = \\begin{pmatrix}\\mu\_1 \\\\ \\mu\_2 \\end{pmatrix} \\in \\mathbb{R}^2](https://latex.codecogs.com/png.latex?%5Cmu%20%3D%20%5Cbegin%7Bpmatrix%7D%5Cmu_1%20%5C%5C%20%5Cmu_2%20%5Cend%7Bpmatrix%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E2 "\mu = \begin{pmatrix}\mu_1 \\ \mu_2 \end{pmatrix} \in \mathbb{R}^2")
and covariane matrix
![\\mathbb{R}^{2 \\times 2} \\ni \\Sigma = \\begin{pmatrix} \\sigma\_{11} & \\sigma\_{12} \\\\ \\sigma\_{2,} & \\sigma\_{22} \\end{pmatrix} \\succ 0](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E%7B2%20%5Ctimes%202%7D%20%5Cni%20%5CSigma%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Csigma_%7B11%7D%20%26%20%5Csigma_%7B12%7D%20%5C%5C%20%5Csigma_%7B2%2C%7D%20%26%20%5Csigma_%7B22%7D%20%5Cend%7Bpmatrix%7D%20%5Csucc%200 "\mathbb{R}^{2 \times 2} \ni \Sigma = \begin{pmatrix} \sigma_{11} & \sigma_{12} \\ \sigma_{2,} & \sigma_{22} \end{pmatrix} \succ 0")
has the following density:

![p(x \| \\mu, \\Sigma) = \\cfrac{1}{2\\pi\|\\Sigma\|} \\exp\\left( -\\cfrac{1}{2}\\left(\\begin{pmatrix}x\_1 \\\\ x\_2 \\end{pmatrix} - \\begin{pmatrix}\\mu\_1 \\\\ \\mu\_2 \\end{pmatrix} \\right)^T \\begin{pmatrix} \\sigma\_{11} & \\sigma\_{12} \\\\ \\sigma\_{21} & \\sigma\_{22} \\end{pmatrix} \\left(\\begin{pmatrix}x\_1 \\\\ x\_2 \\end{pmatrix} - \\begin{pmatrix}\\mu\_1 \\\\ \\mu\_2 \\end{pmatrix} \\right) \\right)](https://latex.codecogs.com/png.latex?p%28x%20%7C%20%5Cmu%2C%20%5CSigma%29%20%3D%20%5Ccfrac%7B1%7D%7B2%5Cpi%7C%5CSigma%7C%7D%20%5Cexp%5Cleft%28%20-%5Ccfrac%7B1%7D%7B2%7D%5Cleft%28%5Cbegin%7Bpmatrix%7Dx_1%20%5C%5C%20x_2%20%5Cend%7Bpmatrix%7D%20-%20%5Cbegin%7Bpmatrix%7D%5Cmu_1%20%5C%5C%20%5Cmu_2%20%5Cend%7Bpmatrix%7D%20%5Cright%29%5ET%20%5Cbegin%7Bpmatrix%7D%20%5Csigma_%7B11%7D%20%26%20%5Csigma_%7B12%7D%20%5C%5C%20%5Csigma_%7B21%7D%20%26%20%5Csigma_%7B22%7D%20%5Cend%7Bpmatrix%7D%20%5Cleft%28%5Cbegin%7Bpmatrix%7Dx_1%20%5C%5C%20x_2%20%5Cend%7Bpmatrix%7D%20-%20%5Cbegin%7Bpmatrix%7D%5Cmu_1%20%5C%5C%20%5Cmu_2%20%5Cend%7Bpmatrix%7D%20%5Cright%29%20%5Cright%29 "p(x | \mu, \Sigma) = \cfrac{1}{2\pi|\Sigma|} \exp\left( -\cfrac{1}{2}\left(\begin{pmatrix}x_1 \\ x_2 \end{pmatrix} - \begin{pmatrix}\mu_1 \\ \mu_2 \end{pmatrix} \right)^T \begin{pmatrix} \sigma_{11} & \sigma_{12} \\ \sigma_{21} & \sigma_{22} \end{pmatrix} \left(\begin{pmatrix}x_1 \\ x_2 \end{pmatrix} - \begin{pmatrix}\mu_1 \\ \mu_2 \end{pmatrix} \right) \right)")

There are no functions in the standard R library to sample from a
bivariate normal. It’s however possible to sample the bivariate normal
using univarite normal distributions through a Gibbs Sampler. It is
because we can show that:

For example, we can use the below Gibbs Sampler implementaiton to sample
from
![\\text{N}\\left( \\begin{pmatrix}1 \\\\ 1 \\end{pmatrix}, \\begin{pmatrix} 1 & 0.5 \\\\ 0.5 & 1 \\end{pmatrix} \\right)](https://latex.codecogs.com/png.latex?%5Ctext%7BN%7D%5Cleft%28%20%5Cbegin%7Bpmatrix%7D1%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%2C%20%5Cbegin%7Bpmatrix%7D%201%20%26%200.5%20%5C%5C%200.5%20%26%201%20%5Cend%7Bpmatrix%7D%20%5Cright%29 "\text{N}\left( \begin{pmatrix}1 \\ 1 \end{pmatrix}, \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix} \right)")

``` r
# Gibbs Sampler for  Example 4

# sample x1
sample_x1 = function() {
  m1 = mu[1] + Sigma[1,2] * (x2 - mu[2]) / Sigma[2,2]
  s1.2 = (Sigma[1, 1] * Sigma[2, 2] - Sigma[1,2]^2) / Sigma[2,2]
        
  x1 <<- rnorm(1, m1, sqrt(s1.2))
}

# sample x2
sample_x2 = function() {
  m2 = mu[2] + Sigma[1,2] * (x1 - mu[1]) / Sigma[1,1]
  s2.2 = (Sigma[1, 1] * Sigma[2, 2] - Sigma[1,2]^2) / Sigma[1,1]
        
  x2 <<- rnorm(1, m2, sqrt(s2.2))
}

# one Gibbs iteration consisits of sample_x1() and sample_x2
gibbs_iter = function() {
  sample_x1()
  sample_x2()
}
```

``` r
# Define mean and covariance matrix of bivariate normal distribution
Sigma <- matrix(1, nrow = 2, ncol = 2) 
Sigma[1,1] = 1; Sigma[1,2] = 0.5; Sigma[2,1] = 0.5; Sigma[2,2] = 1
mu = c(1, 1)

# Run Gibbs until mixing
x1 = -1; x2 = 1 # set initial point
niter = 1000
for(i in 1:niter) {
  gibbs_iter()
}
```

``` r
N = 1000
sample = matrix(0, N, 2)

# get bivariate normal sample from the converged gibbs
for(n in 1:1000) {
  gibbs_iter()
  sample[n,] = c(x1, x2)
}

# plot out the data
ggplot(data=data.frame(x1 = sample[,1], x2 = sample[,2])) + geom_point(aes(x1,x2)) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Detour: OOP in R

Reference classes, called RC for short, are the newest OO system in R.
They were introduced in version 2.12. RC implements message-passing OO,
so methods belong to classes, not functions. $ is used to separate
objects and methods, so method calls look like person$say(“hello”). RC
objects are also mutable: they don’t use R’s usual copy-on-modify
semantics, but are modified in place.

We you use setRefClass() to create a new RC class. A class has its
fields (attributes) and methods. In addition, RC uses the keyword
“contains” to indicate the parent class from which the current class
inherits from.

``` r
Person <- setRefClass("Person",
    # fields
    fields = list(name = "character", age = "numeric",gender= "character"),
    #  methods
    methods = list(
      # Constructors
      initialize = function(name_init = "", age_init = 0, gender_init = "") {
        name <<- name_init
        age <<- age_init
        gender <<- gender_init
      },
      greet = function() {
        print(paste0("Hello, my name is ", name))
      }
    )
)
```

``` r
Student <- setRefClass("Student",
    # use keyword contains to inherit from the parent class
    contains = "Person",
    fields = list(gpa="numeric"),
    methods = list(
      # Constructor
      initialize = function(name_init = "", age_init = 0, gender_init = "", gpa_init = 0) {
        name <<- name_init
        age <<- age_init
        gender <<- gender_init
        gpa <<- gpa_init
      },
      report_gpa = function() {
        print(paste0("My gpa is ", gpa))
      }
    )
)
```

``` r
john = Student$new("John", 21, "M", 3.50) # intialize an instance of the Student class
john$report_gpa() # call method report_gpa()defined in Student
```

    ## [1] "My gpa is 3.5"

``` r
john$greet() # call greet()inherited from the parent Person
```

    ## [1] "Hello, my name is John"

Note that RC objects are mutable, i.e., they have reference semantics,
and are not copied-on-modify.

``` r
mary = john
mary$name = "Mary"
print(john$name)
```

    ## [1] "Mary"

As a consequence, you should use copy() method to clone an object.

``` r
john$name = "John"
mary = john$copy()
mary$name = "Mary"
mary$gender = "F"
print(john$name)
```

    ## [1] "John"

``` r
print(john$gender)
```

    ## [1] "M"

``` r
print(mary$name)
```

    ## [1] "Mary"

``` r
print(mary$gender)
```

    ## [1] "F"

## Example 5

In the last example, we’ll demonstrate how to use Gibbs Sampler to
sample the posterior for data coming from a mixture of two normal
distribution
![p(x) = 0.5 \\cdot \\text{N}(x \| -5, 1) + 0.5 \\cdot \\text{N}(x \| 3, 1)](https://latex.codecogs.com/png.latex?p%28x%29%20%3D%200.5%20%5Ccdot%20%5Ctext%7BN%7D%28x%20%7C%20-5%2C%201%29%20%2B%200.5%20%5Ccdot%20%5Ctext%7BN%7D%28x%20%7C%203%2C%201%29 "p(x) = 0.5 \cdot \text{N}(x | -5, 1) + 0.5 \cdot \text{N}(x | 3, 1)")
(<span style="color:blue">dist 1</span>).

<span style="color:red">**Exercise**</span> : show that <span
style="color:blue">dist 1</span> is equivalent to <span
style="color:blue">dist 2</span> defined as:

The model equivalence justifies the following codes to generate data:

``` r
n = 1000
x = rep(0, n) # initialize a sample of size n
mu = c(-5, 3)
  
for(i in 1:n) {
  zi = rbinom(1, 1, 0.5) + 1 # z[i] specifies which cluster x[i] belongs to
  x[i] = rnorm(1, mu[zi], 1)
}
  
ggplot(data.frame(x = x), aes(x=x)) + geom_histogram(color="black", fill="white", binwidth = 0.25) + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Given a data sample and assume that we only know variances, define the
model as follows:

Under this model, we can show that

where

``` r
Normal2DGibbs <- setRefClass("Normal2DGibbs",
    fields = list(mu = "numeric", pi = "numeric", z = "numeric", sigma.2 = "numeric", n = "numeric", x = "numeric"),
    methods = list(
      initialize = function(mu1_init, mu2_init, pi_init, s.2, x) {
        mu <<- c(mu1_init, mu2_init)
        pi <<- pi_init
        sigma.2 <<- s.2
        
        x <<- x
        n <<- length(x)
      },
      sample_z = function() {
        for(i in 1:n) {
          p0 = (1 - pi) * dnorm(x[i], mu[1], sqrt(sigma.2))
          p1 = pi * dnorm(x[i], mu[2], sqrt(sigma.2)) 
          
          z[i] <<- rbinom(1, 1, p1 / (p0 + p1))
        }
      },
      sample_mu = function() {
        for(k in 0:1) {
          nk = sum(z == k)
          tauk =  1 / (1 + nk)
          mk = tauk * sum(x[z == k])
          
          mu[k + 1] <<- rnorm(1, mk, sqrt(tauk))
        }
        
        mu <<- mu[order(mu)]
      },
      sample_pi = function() {
        pi <<- rbeta(1, 1 + sum(z == 0), 1 + sum(z == 1))
      },
      iter = function() {
        sample_z()
        sample_pi()
        sample_mu()
      },
      save_to_disk = function(path = 'Mix2NormGibbs.rds') {
        obj = .self
        saveRDS(obj, file = path)
      }
    )
)
```

``` r
gibbs = Normal2DGibbs(-1, 1, 0.1, 1, x) # intialize gibbs

niters = 1000
for(i in 1:niters) {
  gibbs$iter() # run gibbs until convergence
}

gibbs$save_to_disk() # save converged gibbs to hard disk
```

``` r
# load converged Gibbs from hard drive
gibbs = readRDS('Mix2NormGibbs.rds')

# N = sample size
N = 1000
mu1 = rep(0, N)
mu2 = rep(0, N)
pi = rep(0, N)

# sample from converged gibbs
for(n in 1:N) {
  gibbs$iter()
  mu1[n] = gibbs$mu[1]
  mu2[n] = gibbs$mu[2]
  pi[n] = gibbs$pi
}
```

``` r
# plot the posterior
ggplot(data.frame(mu1 = mu1), aes(x= mu1)) + geom_histogram(color="black", fill="white", binwidth = 0.01) + ggtitle("Posterior mu1 | x") +  theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
ggplot(data.frame(mu2 = mu2), aes(x= mu2)) + geom_histogram(color="black", fill="white", binwidth = 0.01) + ggtitle("Posterior mu2 | x") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

``` r
ggplot(data.frame(pi = pi), aes(x= mu2)) + geom_histogram(color="black", fill="white", binwidth = 0.025) + ggtitle("Posterior pi | x") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->

# MCMC Diagnostic

In the last the two examples, how do we ensure our Gibbs Sampler
actually converge ? In theory, MCMC are guaranteed to converge to its
stationary distribution as iterations go to infinity. We can not of
course run our markov chain forever. While there are no ways to tell
with absolute certainty that the Gibbs Sampler has finally converged,
there are tools to give us some assurnaces. In this section, we’ll see
how to use traceplots, acf plots and effective sample size for MCMC
diagnosis.

## Multiple chains & Trace plots

In practice, it’s critical to run multiple Gibbs Samplers and cross
compare their mixing results. Let’s redo Example 4 but run four separate
Gibbs chains this time. If all of the chains converge, they should cover
the same range of values for each of the model parameters.

``` r
# create 4 Gibbs Sampler chains
chain1 = Normal2DGibbs(-1, 1, 0.1, 1, x)
chain2 = Normal2DGibbs(-3, 2, 0.1, 1, x)
chain3 = Normal2DGibbs(-10, 10, 0.1, 1, x)
chain4 = Normal2DGibbs(-1, 5, 0.1, 1, x)

# mu & pi tracking time series
mu1_chain1 = rep(0, N)
mu1_chain2 = rep(0, N)
mu1_chain3 = rep(0, N)
mu1_chain4 = rep(0, N)
mu2_chain1 = rep(0, N)
mu2_chain2 = rep(0, N)
mu2_chain3 = rep(0, N)
mu2_chain4 = rep(0, N)
pi_chain1 = rep(0, N)
pi_chain2 = rep(0, N)
pi_chain3 = rep(0, N)
pi_chain4 = rep(0, N)

# run each chain for niters
niters = 1000
for(i in 1:niters) {
  chain1$iter()
  mu1_chain1[i] = chain1$mu[1]
  mu2_chain1[i] = chain1$mu[2]
  pi_chain1[i] = chain1$pi
  
  chain2$iter()
  mu1_chain2[i] = chain2$mu[1]
  mu2_chain2[i] = chain2$mu[2]
  pi_chain2[i] = chain2$pi
  
  chain3$iter()
  mu1_chain3[i] = chain3$mu[1]
  mu2_chain3[i] = chain3$mu[2]
  pi_chain3[i] = chain3$pi
  
  chain4$iter()
  mu1_chain4[i] = chain4$mu[1]
  mu2_chain4[i] = chain4$mu[2]
  pi_chain4[i] = chain4$pi
}
```

``` r
# trace plot for mu1
ggplot() + 
  geom_line(data = data.frame(x = 1:niters, y = mu1_chain1), aes(x,y),
                     color = "red") + 
  geom_line(data = data.frame(x = 1:niters, y = mu1_chain2), aes(x,y),
                     color = "blue") +
  geom_line(data = data.frame(x = 1:niters, y = mu1_chain3), aes(x,y),
                     color = "green") +
  geom_line(data = data.frame(x = 1:niters, y = mu1_chain4), aes(x,y),
                     color = "orange") +
  xlab("Iterations") + ylab("mu1") + ggtitle("Trace plot for mu1") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
# trace plot for mu1
ggplot() + 
  geom_line(data = data.frame(x = 1:niters, y = mu2_chain1), aes(x,y),
                     color = "red") + 
  geom_line(data = data.frame(x = 1:niters, y = mu2_chain2), aes(x,y),
                     color = "blue") +
  geom_line(data = data.frame(x = 1:niters, y = mu2_chain3), aes(x,y),
                     color = "green") +
  geom_line(data = data.frame(x = 1:niters, y = mu2_chain4), aes(x,y),
                     color = "orange") +
  xlab("Iterations") + ylab("mu2") + ggtitle("Trace plot for mu2") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
# trace plot for pi
ggplot() + 
  geom_line(data = data.frame(x = 1:niters, y = pi_chain1), aes(x,y),
                     color = "red") + 
  geom_line(data = data.frame(x = 1:niters, y = pi_chain2), aes(x,y),
                     color = "blue") +
  geom_line(data = data.frame(x = 1:niters, y = pi_chain3), aes(x,y),
                     color = "green") +
  geom_line(data = data.frame(x = 1:niters, y = pi_chain4), aes(x,y),
                     color = "orange") +
  xlab("Iterations") + ylab("pi") + ggtitle("Trace plot for pi") + theme_bw()
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

## Sample autocorrelation

In theory, after the markov chain converges, their sample will be then
iid. Therefore, one way to check whether these chains have mixed is to
draw sample and either use acf() or compute effective sample size.

``` r
# Obtain a sample from each chain
# N = sample size
N = 1000

mu1 = rep(0, N)
mu2 = rep(0, N)
pi = rep(0, N)

# draw sample from all four converged chains
for(n in 0: (N/4 - 1))  {
  chain1$iter()
  mu1[4*n + 1] = chain1$mu[1]
  mu2[4*n + 1] = chain1$mu[2]
  pi[4*n + 1] = chain1$pi
  chain2$iter()
  mu1[4*n + 2] = chain2$mu[1]
  mu2[4*n + 2] = chain2$mu[2]
  pi[4*n + 2] = chain2$pi
  chain3$iter()
  mu1[4*n + 3] = chain3$mu[1]
  mu2[4*n + 3] = chain3$mu[2]
  pi[4*n + 3] = chain3$pi
  chain4$iter()
  mu1[4*n + 4] = chain3$mu[1]
  mu2[4*n + 4] = chain3$mu[2]
  pi[4*n + 4] = chain3$pi
}
```

Ideally, for acf plot, at
![\\alpha = 5](https://latex.codecogs.com/png.latex?%5Calpha%20%3D%205 "\alpha = 5")%
significant level, we expect only 1/20 lag to be outside the confidence
interval (blue lines). The acf plot for chains’ sample look good enough.

``` r
require(forecast)
# plot out the acf plot
ggAcf(mu1)
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
ggAcf(mu2)
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

``` r
ggAcf(pi)
```

![](551_gibbs_files/figure-gfm/unnamed-chunk-33-3.png)<!-- -->

As mentioned eariler, we can also check the effective sample size. If
chains have not yet mixed, their sample will be highly correlated and
their effective sample size (ESS) will be small. In our case, the ESS
out of 1000 Gibbs iterations are sufficiently large for all the
parameters .

``` r
require(coda)
print(paste0("mu1 ESS is ", effectiveSize(mu1)))
```

    ## [1] "mu1 ESS is 767.068076612478"

``` r
print(paste0("mu2 ESS is ", effectiveSize(mu2)))
```

    ## [1] "mu2 ESS is 797.686198384314"

``` r
print(paste0("pi ESS is ", effectiveSize(pi)))
```

    ## [1] "pi ESS is 642.547648397276"

In summary, the traceplots, acf and effective size diagnosis gave us
evidences that our Gibbs chains have converged.

**Reference:**

\[1\] Long Nguyen, Stats 551 lecture notes and class material, 2021.

\[2\] Peter D. Hoff, A First Course in Bayesian Statistical Methods.
Springer, 1st ed, 2009.

\[3\] Hadley Wickham, Advanced R. Chapman & Hall, 2nd ed, 2019.
