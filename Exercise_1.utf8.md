---
title: "Project 1 - SPATIAL STATISTICS"
author: "Kwaku Peprah Adjei, Amir Ahmed"
#date: "1/23/2020"
output:
  pdf_document:
    fig_caption: yes
    includes:
      in_header: preamble.tex
  html_document:
    df_print: paged
    includes:
      in_header: preamble.tex
editor_options:
  chunk_output_type: console
---





\newpage
# Problem 1
Let $\lbrace r(x) : x \in \text{D} : \left[1, 50\right]\subset\mathbb{R}^1\rbrace$. Assume it is modelled as a stationary one dimensional Gaussian random field with the following parameters:
\begin{equation}\label{prob1eq}
    \begin{split}
    E\lbrace r(x) \rbrace &= \mu_r = 0 \\
    Var\lbrace r(x)\rbrace &= \sigma_r^2 \\
    Corr\lbrace r(x), r(x')\rbrace &= \rho_r(\tau) \\
    \tau = \dfrac{|x-x'|}{10}
    \end{split}
\end{equation}
where $\rho_r(\tau)$ is the spatial correlation function. Discretize $D$ as $L = \lbrace 1, 2, \dots, 50 \rbrace$, we will later look at the discretized Gaussian random field $\lbrace r(x) : x\in L \rbrace$.


## Problem 1a
To ensure the positive definetness of all covariance matrices, we rely on having a positve definite corrolation function. A positive definite correlation function guarantees that the covariance matrix of our Gaussian random field is well defined for all choices of observation points and for all dimensions.  

A function $c(\vect \tau): \mathbb{R}^q$ i positive definite if the following is satisfied:
\begin{equation}\label{posdefeq}
    \begin{split}
        \sum_{i=1}^n\sum_{j=1}^n\alpha_i\alpha_j c(\vect x_i - \vect x_j) \geq 0,
    \end{split}
    \quad
    \begin{split}
        &\forall \alpha_1, \dots, \alpha_n \in \mathbb{R} \\
        &\forall n \in \mathbb{N}, \quad n \geq 2 \\
        &\forall \vect x_1, \dots \vect x_n \in \mathbb{R}^q
  \end{split}
\end{equation}
Any $n \times n$ matrix, $\matr Q = \left[ c(\vect x_i - \vect x_j) \right]_{ij}$ constructed with a positive definite $c(\vect \tau)$ with vectors $\lbrace \vect x_i \rbrace_{i=1}^n$ chosen as in \eqref{posdefeq} would clearly be positive definite. As for all $\vect \alpha = (\alpha_1, \dots, \alpha_n) \in \mathbb{R}^n$ we have:
\begin{equation}
    \vect \alpha^T \matr Q \vect \alpha =   \sum_{i=1}^n\sum_{j=1}^n\alpha_i\alpha_j c(\vect x_i - \vect x_j) \geq 0
\end{equation}
So it would be safe to use $\matr Q$ as a covariance matrix, and it would satisfy the requirements of the multivariate normal distribution.

We now go on to look at the the powered exponential correlation function:

\begin{equation}
    \begin{split}
            \rho(\tau) = \exp(-\tau^{\nu}) 
    \end{split}, \quad
    \begin{split}
    \nu \in \left(0, 2 \right]
    \end{split}
\end{equation}
for later case studies we will look at parameters $\nu_r \in \lbrace 1, 1.9 \rbrace$. We note that and increased $\nu$ ($\nu >> 1$) would lower the exponent and thus give less correlation. 

We also want to study the Matern correlation function
\begin{equation}
    \begin{split}
            \rho_r(\tau) = \dfrac{2^{1-\nu}}{\Gamma(\nu)}\tau^v\mathcal{B}_\nu(\tau), \quad
    \end{split}
    \begin{split}
        \nu \in \mathbb{R}_+
    \end{split}
\end{equation}
$\mathcal{B}_\nu$ is the bessel function. We want to look at the matern with parameters $\nu_r \in \lbrace 1, 3 \rbrace$. 

For both correlation functions use $\sigma^2_r \in \lbrace 1, 5\rbrace$

![\label{fig:fig1a1} Plotting matern and exponential correlation](Exercise_1_files/figure-latex/fig1a1-1.pdf) 


The two correlation functions are plotted in Figure \ref{fig:fig1a1} for $\tau \in \mathbb{R}_\oplus$. We note that these functions both are positive definite. From the figures we also note that both correlation functions for the different parameters seems to satisfy ergodicity, correlation drops to zero the further apart to points are, this is a neccesary trait in our Gaussian RF. 

For the matern corrolation we see that an increase of $\nu$ lead to slower fall in correlation with distance. We also note that the Matern seems to drop slower in correlation than the exponential correlation function for the given parameters.

For higher values of $\nu$ and  when $\tau >> 1$ the powered exponential seems to drop faster in correlation than for higher $\nu$. The oppsite seems to be the trend when $\tau << 1$

The variogram function is defined as:
\begin{equation}\label{eq:variogram}
    \begin{split}
        \gamma_r(\tau)  &= \frac{1}{2} Var\lbrace r(x) - r(x') \rbrace \\
        &= \frac{1}{2} (Var\lbrace r(x) \rbrace + Var\lbrace r(x') \rbrace - 2 Corr\lbrace r(x), r(x')) \rbrace\\
        &=\sigma_r^2(1-\rho_r(\tau))
    \end{split}
\end{equation}

Looking at \eqref{eq:variogram} we see that if $\sigma_r^2 = 1$ then the corrolation function and the variogram functions would respectively increase and decrease at same rate in with respect to an increase in $\tau$. In essence the variogram function tells us how much variance we have in our estimation of $x'$ when it is $\tau$ away from a observed point $x$

![\label{fig:fig1a2} Plotting variograms](Exercise_1_files/figure-latex/fig1a2-1.pdf) 

We display the variograms for our model parameters in Figure \ref{fig:fig1a2}.

In both figures we see that the value of $\sigma^2$ sets a roof on the variance of points far away. Higher values of $\sigma^2$ also seems to increase the max achieved variance. For the matern variogram increased $\nu$ seems to increase corrolation to neighbouring points. 


\newpage

## Problem 1b                    

We use the corrolation function $\rho_r(\tau)$ to construct a corrolation matrix. The corrolation matrix would have the following form:
\begin{equation}
    \matr \Sigma_r^\rho = 
    \begin{bmatrix}
        \rho_r(1, 1) & \rho_r(1, 2) & \dots & \rho_r(1, 50) \\
        \rho_r(2, 1) & \rho_r(2, 2) & \dots & \rho_r(2, 50) \\
        \vdots & \vdots & \ddots & \vdots \\
        \rho_r(50, 1) & \rho_r(50, 2) & \dots & \rho_r(50, 50)
    \end{bmatrix}
\end{equation}
$\rho_r(x, x')$ denotes $\rho_r(\tau) = \rho_r(|x-x'|/10)$. With variance $\sigma_r^2$ the prior get the covariance matrix $\matr \Sigma_r = \sigma_r^2 \matr \Sigma_r^\rho$. From \eqref{prob1eq} we have an expected value of:
\begin{equation}
    \vect \mu_r = (0, \dots, 0)^T
\end{equation}
Which gives the prior a distribution of:
\begin{equation}
    \phi_{50}(\vect r ; \vect \mu_r, \matr \Sigma_r)
\end{equation}
Which has pdf:
\begin{equation}
    (2\pi)^{-\frac{50}{2}}\det(\matr \Sigma_r)^{-\frac{1}{2}}\exp\left(-\frac{1}{2}\vect x^T\matr \Sigma_r^{-1}\vect x\right)
\end{equation}

We simulate four realisations of the field for the different parameters.

\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig1b1-1} 

}

\caption{\label{fig:fig1b1} Realizisations of prior, matern}\label{fig:fig1b1}
\end{figure}

\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig1b2-1} 

}

\caption{\label{fig:fig1b2} Realizisations of prior, powered exponential}\label{fig:fig1b2}
\end{figure}


The results are plotted in Figure \ref{fig:fig1b1} and Figure \ref{fig:fig1b2}. We note that an increased variance seems to increase fluctuation (comparing c) d) g) and h) to the others). For the matern we also see that an increased $\nu$ seems to make the process smoother. For both correlation functions increased $\nu$ seems to increase "air-time" when the process is strays away from expected 0. 


\newpage

## Problem 1c

We now want to develop a posterior model. Assume that we have observed the values at $x \in \lbrace 10, 25, 30 \rbrace$. Organise them as:
\begin{equation}
    \lbrace d(x); x \in \lbrace 10, 25, 30 \rbrace \subset L  \rbrace 
\end{equation}
We also assume we have an observation error $\sigma^2_\epsilon = \lbrace 0, 0.25 \rbrace$. In general we have
\begin{gather}
    d(x) = r(x) + \epsilon(x), \quad x \in \lbrace 10, 25, 30 \rbrace \\
   \epsilon(x) \sim N(0, \sigma_\epsilon^2), \iid \\
   r(x) \text{ and } \epsilon(x) \text{ independent} \\
   \epsilon(x) \text{ and } \epsilon(x') \text{ independent identically distributed} 
\end{gather}

As we both $\epsilon(x)$ and $r(x)$ are Gaussian, a linear product of the two would also be. 

We further note:
\begin{equation}
    E(d(x)) = E(r(x)) + E((\epsilon(x))) = 0
\end{equation}
When $x\neq x'$
\begin{equation*}
    \begin{split}
        \Cov(d(x), d(x')) &= \Cov(r(x)+\epsilon(x), r(x') + \epsilon(x')) \\
        &= \Cov(r(x), r(x')) + \Cov(\epsilon(x), \epsilon(x')) \\
        &= \sigma_r^2\rho_r(\tau)
    \end{split}
\end{equation*}
And: 
\begin{equation}
    \begin{split}
        \Cov(d(x), d(x)) &= \Cov(r(x)+\epsilon(x), r(x) + \epsilon(x)) \\ 
        &= \Cov(r(x), r(x)) + \Cov(\epsilon(x), \epsilon(x)) \\
        &= \sigma_r^2\rho_r(\tau) + \sigma_e^2
    \end{split}
\end{equation}
We organise the observed points in a $k \times 1$ vector,$\vect y = (10, 25, 30)^T$, where $k = 3$. Further let the observed values of the random field be: 
\begin{equation}
    \vect d(\vect y) = ( d(10), d(25), d(30))^T
\end{equation}
If we denote $\matr \Sigma_\rho^d$ as the correlation matrix generated by $\rho(\tau)$ between the points $\vect x$
we see:
\begin{equation}
    \Cov(\vect d(\vect x), \vect d(\vect x')) = \sigma^2_e \matr I_k +\sigma^2_r\matr \Sigma^\rho_d
\end{equation}

Using the calculations above $\vect d$ thus have the following pdf:
\begin{equation}
    p(\vect d(\vect x)|\sigma_r^2, \sigma_e^2) = (2\pi)^{k/2}\det(\matr \Sigma_\rho^d)^{-1/2}\exp(\frac{1}{2}\vect d^T(\matr\Sigma_\rho^d)^{-1}\vect d)
\end{equation}
The value $\sigma_e^2$ would then have the following likelihood function:
\begin{equation}
    L(\sigma_e^2 | \vect d(\vect x), \sigma_r^2)=p(\vect d(\vect x)|\sigma_r^2, \sigma_e^2)
\end{equation}
The integral 
\begin{equation}
    \int_{0}^{\infty}L(\sigma_e^2 | \vect d(\vect x), \sigma_r^2)d\sigma^2 = \int_{0}^{\infty}p(\vect d(\vect x)|\sigma_r^2, \sigma_e^2) d\sigma^2
\end{equation}
goes to infinity, thus $L(\sigma_e^2 | \vect d(\vect x), \sigma_r^2)$ is not a pdf. If we try to estimate $\sigma^2_\epsilon$ by the expected value of $\sigma_e^2$ derived when believing $L(\sigma^2_\epsilon|\cdot)$ is a pdf, we would expect values that goes to infinity, which would be wrong. 



## Problem 1d

Assume now that $\sigma_r^2 = 5$ and that we have observed $\vect d(\vect y)$ with error $\sigma_\epsilon^2 = \in \lbrace 0, 0.25 \rbrace$. Want to find the distribution of $\left[ \vect r | \vect d \right]$. Know that both $\vect r$ and $\vect d$ are multivariate normal. Thus conditioning $\vect r$ on $\vect d$ would give a multivariate normal distribution with the following parameters: 

\begin{equation}
     \vect\mu_{\vect l | \vect d} = \vect \mu_r+ \matr \Sigma_{\vect r, \vect d} \matr \Sigma_{\vect d}^{-1}(\vect d - \vect \mu_{\vect d})
\end{equation}

\begin{equation}
      \vect\Sigma_{\vect r | \vect d} = \matr \Sigma_{\vect r} - \matr \Sigma_{\vect r, \vect d} \matr \Sigma_{\vect d}^{-1}\Sigma_{\vect d, \vect r}
\end{equation}

Note $\Sigma_{\vect r, \vect d} = \Sigma_{\vect d, \vect r}^T$ and they are respectively $50\times 3$ and $3\times 50$ where the first has elements $[\Cov\lbrace r(x_i), d(x_j)\rbrace]_{ij}$ where:

\begin{equation}
    \begin{split}
        \Cov\lbrace r(x_i), d(x_j)\rbrace &= \Cov\lbrace r(x_i), r(x_j) + \epsilon(x_j)\rbrace  \\
        &= \Cov\lbrace r(x_i), r(x_j)\rbrace \\ 
        &= \sigma_r^2\rho_r(x_i, x_j)
    \end{split}
\end{equation}

Assuming we have realisations of $\vect d$ from simulated data we can study how $\vect r | \vect d$ acts. 

\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig1d1-1} 

}

\caption{\label{fig:fig1d1} Posterior estimation, matern $\sigma^2 = 5$, $\nu = 3$}\label{fig:fig1d1}
\end{figure}



We compute two corresponding predictions for the spatial variable $\lbrace \hat r(x);x \in L\rbrace$ with associated 0.9 predictions intervals and display them in Figure \ref{fig:fig1d1}. We use the same realisations Trial 1 in d) in Figure \ref{fig:fig1b1}. 

Here we have parameters $\nu_r  = 3$ and  $\sigma_r^2 = 5$. As discussed earlier, we would with these values see much corrolation. From the variogram, points with distance up to about 20 in distance would impact each other, this is reflected in the figure, as the observation points shifts large parts of the process away from 0. The high correlation decreases the width of the 90\% confidence intervals close to the observations. In the interval $x\in (25,30)$ the 90\% band is very narrow in the exact observation case, comapred to the one with observation error. Comparing no variance in observation to variance the points close to the observation points in no observation error is of course more accurate, but further away i.e. at point 20 the increased variance seem to have little impact on the confidence band. The real process behind, keeps within 90\% although we would expect more of them to stray of, as we have 47 non observed-points.


## Problem 1e
We now simulate 100 realisation's with the same parameters as in 1d) and calculate empirical mean and empirical 90\% confidence bands. 
\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig1e1-1} 

}

\caption{\label{fig:fig1e2} Posterior estimation 100 realisation, matern $\sigma^2 = 5$, $\nu = 3$}\label{fig:fig1e1}
\end{figure}

The results are illustrated in Figure \ref{fig:fig1e2}. We also include the theoretical plots that the simulations are drawn from. For both with observation error and without observations error the empirical estimations seems to do well, and overlaps quite well with the theoretical bands. The largest differences between the two are seen furthest away from the observations points. Between $x\in(10,25)$ theoretical bands seems to have a slight tendency to be wider that the empirical bands, however the trend is not too clear. Close to the observation points the bands seems to be equally wide without observation error, with obserbation error there seems to be slightly more variation.  



## Problem 1f

Want to evaluate 

\begin{equation}
    A_r = \int_D I(r(x) > 2)dx
\end{equation}

\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig1f1-1} 

}

\caption{\label{fig:fig1f1} Posterior estimation 100 realisation, matern $\sigma^2 = 5$, $\nu = 3$}\label{fig:fig1f1}
\end{figure}

The situation is illustrated in Figure \ref{fig:fig1f1}


A prediction of Â using our 100 realizations from the posterior model with $\sigma_\epsilon^2$ is:
\begin{equation}
    \hat A_r = \frac{1}{100}\sum_{i=1}^{100}\sum_{x\in L}I(r_i(x)>2)
\end{equation}
Applying this to our realisation we get: Â = 30.16 and Var(Â ) = 32.519596

Where we use the approximation:
\begin{equation}
    A_r = \int_D I(r(x) > 2)dx \approx \sum_{x \in L}I(r(x)>2) dx 
\end{equation}
Another estimator is: 
\begin{equation}
    \tilde A_r = \sum_{x \in L} I(\hat r(x) > 2)
\end{equation}
With our data we get: Ã =  29

We note that in our case $Â_r > Ã_r$. 


If we let:
\begin{equation}
    g(r(x)) = \int_D I(r(x) > 2) dx
\end{equation}
Then $g(x)$ is convex as it is an integral of a positive value. 

Using Jensen's Inequality we thus get: 
\begin{equation}
 Eg(r(x)) \geq g(Er(x)) = \int_D I(Er(x) > 2) =\int_D I(\hat r(x) > 2)
\end{equation}
The same would be the case for the discretized stimators, i.e. we have:
$$E\hat A_r \geq \tilde A_r$$
And is why we expect $\hat A_r > \tilde A_r$. Which is also what we observed. 

Another way to look at the problem is the following:

Let $Y_{ij}$ be that event that the j-th observation at the i-th measurement point ($0 \leq i \leq 50$) is over or under 2 in the posterior model. We define:
\begin{equation}
    \begin{split}
        & Y_{ij} := I(r(X_{ij})> 2)
    \end{split}
\end{equation}
Where $X_{ij}$ is the j-th observation positionat the i-th realisation of the posterior. As the 100 posterior realisations are identically distributed  and independent. Each $Y_{ij}$ is independent grouped for position j.  Then $Y_{ij}$ is Bernoulli distributed. With probability $p_{i}$ of being higher than 2. Where:

\begin{equation}
    p_i = 1 - \Phi\left(\frac{2 -  \mu_i}{\sigma_i}\right)
\end{equation}
Where $\mu_i = \left[\vect \mu_{l|d}\right]_i$ and $\sigma_i^2 = \left[\matr \Sigma_{l|d}\right]_{ii}$

We know that $\hat p_i = 1/100\sum_{i=1}^{100}Y_{ij}$ is an unbiased estimator for $p_i$

The theoretical and estimated values of $p_i$ are illustrated in Figure \ref{referenceHere}, and seem to match the situation and each other quite well. They are scaled such that the height of the plot (-5, 5) would be $p_i = 1$ 

This thinking can also be applied to estimating $A_r$

One alternative $\hat A_r$ is: 
\begin{equation}
    \begin{split}
        \hat A_r =  \frac{1}{100}\sum_i^{50}w\sum_j^{100}I(X_{ij} > 2)  
    \end{split}
\end{equation}
We note that:
\begin{equation}
    \begin{split}
        E\hat A_r &=  \frac{1}{100}\sum_i^{50}w\sum_j^{100}EI(X_{ij} > 2) \\ &=   \frac{1}{100}\sum_{i=1}^{50}w\sum_{j=1}^{100}\left[1 - \Phi\left(\frac{2 -  \mu_i}{\sigma_i}\right)\right] \\
        &= \sum_i^{50}wp_i
    \end{split}
\end{equation}
Whre $w$ is the width between each observation, in our case $w=1$. We see that when number of observations goes to infinity, the above sum goes toward the integral we want to estimate. 

This can also be used to show that we expect $\hat A > \tilde A$. By using that:
\begin{equation}
    P(\hat r(x_i) > 2) = 1 - \Phi\left(\frac{2-\mu_i}{\sigma_i}\sqrt{n}\right) < 1 - \Phi\left(\frac{2-\mu_i}{\sigma_i}\right)
\end{equation}
The inequality can be easily proven by using this. 

# Problem 1g)
We have studied the matern and the powered exponential correlation function, and applied them to simulated data, to estiamte posteriors and their confidence. We find it interesteing how a change in the $\nu$ seem to change the differentiability of the end smoothness. Overall the empirical estimation methods seems to work well and are consistent with theory. 


\newpage
# Problem 2 Gaussian RF - real data
Given the domain $D = [(0,315) \times (0,315)] \subset \mathbb{R^2}$. We let $\vect{d}= r(\vect{x_1^0)}, ..., r(\vect{x_{52}^0)})^T$. 

## Problem 2a: Display of the data

The data is displayed with an image plot, a contour plot and the exact points as shown in the figure \ref{fig:fig1} below.
![\label{fig:fig1}, fig.height=0.5](Exercise_1_files/figure-latex/fig1-1.pdf) 

It was observed that the data points did not move in the same direction as with the x and y cordinates (see figure \ref{fig:fig2} a;b). This suggest that the data is not mean stationary. Moreover, a density plot for the data points in figure \ref{fig:fig2} c) shows a skewness in the data, making the Guassianity assumption doubtful. Hence a stationary Gaussian RF may not be appropriate.

![\label{fig:fig2} Plot of the data points with respect to their a) x-cordinates and b) y-cordinates; and c) shows the density distribution of the terrain elevation observations.](Exercise_1_files/figure-latex/fig2-1.pdf) 


## Problem 2b
Let 
$$
\{r (\vect{x}); \vect{x} \in D \subset \mathbb{R^2}\}
$$
be the Gaussian RF that is used to model the domain $D$.

Given that $E\{r(x)\} = (\vect{gx})^T \vect{\beta_r}$, $Var\{r(\vect{x})\} = \sigma_r^2$ and $Corr(r(\vect{x}), r(\vect{x'})) = \rho_r(\frac{\tau}{\xi})$. We assume that $\sigma_r^2, \xi$ are assumed known but $\vect{\beta_r}$ is unknown. This is therefore a universal krigging problem. 

Let the krigging predictor be:
$$
\vect{\hat{r}_0} = \vect{\alpha}^T \vect{r^d} 
$$

We discretise the predictor as:
$$
\{\vect{r_{\Delta} (x)} = r(\vect{x}) - \mu_r^0 - \sum_{i=1}^{n_g} \beta_r^j g_j(\vect{x}); \vect{x} \in  D\}
$$

For the estimator to be unbiased, 

$$E\{\vect{\hat{r}_0} - \vect{{r}_0} =0 \} \implies \sum_{i=1}^m \alpha_iE\{r_i^d\} - E\{ \vect{{r}_0}\} = 0$$ 

$$\sum_{i=1}^m \sum_{j=1}^{n_g} \alpha_i \beta_r^j g_j(\vect{x}_i) = \sum_{j}\beta_r^j g_j(\vect{x}_0)$$


For the estimator to be unbiased, 

$\sum_{i=1}^m \alpha_i g_j(\vect{x}_i) = \sum_{j}\beta_r^j g_j(\vect{x}_0)$.

\begin{equation*}
    \begin{split}
      Var\{\vect{\hat{r}_0} - \vect{{r}_0} \} &= E\{(\vect{\hat{r}_0} - \vect{{r}_0})^2 \} \\
                                       &= Var\{\alpha_i \{r_i^d\} - \vect{{r}_0}\} \\
                                        &= \sigma^2 \sum_{i=1}^n\sum_{j=1}^m \alpha_i\alpha_j\rho_{ij} + \sigma^2 + 2 \sigma^2\sum_{j=1}^m\alpha_j\rho_{j0}\\
    \end{split}
\end{equation*}


Hence, we find $\vect{\hat\alpha}$ such that
$$
\vect{\hat\alpha} = argmin_{\vect{\alpha}} Var\{\vect{\hat{r}_0} - \vect{{r}_0} \}
$$
and subject to the constraint $\sum_{i=1}^m \alpha_i g_j(\vect{x}_i) = \sum_{j}\beta_r^j g_j(\vect{x}_0)$ for $j= 1,2,...,n_g$.

## Problem 2c

Considering the case with $E(r(\vect x)) = \beta_1$, we estimated the universal krigging predictor and variance as follows:



![\label{fig:fig3} Krigging predictions and prediction variance of the ordinary krigging method.](Exercise_1_files/figure-latex/fig3-1.pdf) 

## Problem 2d

The resulting polynomial function becomes:
$$
\vect{(gx)} = (1, x_v, x_h, x_vx_h, x_v^2, x_h^2)
$$
The expected value of $r(\vect{x})$ then is:
$$
E \{r(\vect{x}) \} = \beta_1 + \beta_2x_v + \beta_3x_h + \beta_4x_vx_h + \beta_5x_v^2 + \beta_6x_h^2 .
$$



![\label{fig:fig4} Predictions and associated variance, fig.height=0.5](Exercise_1_files/figure-latex/fig4-1.pdf) 
We present the predictions and the associated variance in Figure \ref{fig:fig4}

## Problem 2e




We now consider the grid node, $\vect{x_0} = (100,100)$. Using the ordinary krigging, we estimated the predicted mean as 838.6781414 and the predicted variance as 10.044233. Assuming Gaussianity of the data,

\begin{equation*}
    \begin{split}
P(r\{\vect{x_0}\} > 850) &=P\bigg( \frac{r\{\vect{x_0}\}- E(r\{\vect{x_0}\})}{\sqrt{Var(r\{\vect{x_0}\})}}> \frac{850- E(r\{\vect{x_0}\})}{\sqrt{Var(r\{\vect{x_0}\})}} \bigg)\\
&= 1 - \Phi\bigg(\frac{850- E(r\{\vect{x_0}\})}{\sqrt{Var(r\{\vect{x_0}\})}} \bigg)
  \end{split}
\end{equation*}

The resulting probability is 0.13.

To obtain the elevation for which it is 0.90 probability that the true elevation is below it, we used the formular,

\begin{equation*}
    \begin{split}    
P(r\{\vect{x_0}\} >r\{\vect{x_{new}}\} ) &= 0.90\\
r\{\vect{x_{new}}\} &= E(r\{\vect{x_0}\}) + \phi(0.90) \sqrt{Var(r\{\vect{x_0}\})}
  \end{split}
\end{equation*}

We obatained 851.55m to be that elevation that satifies the preamble.



\newpage 

# Problem 3

We consider the stationary Gaussian RF $\lbrace r((\vect x); \vect x \in D \subset \mathbb{R}^2$. With $D:\left[(1,30) \times (1,30) \right]$. 

With: 
\begin{equation}
    \begin{split}
        E\lbrace r(\vect x) \rbrace &= \mu_r = 0 \\
        Var\lbrace r(\vect x) \rbrace &= \sigma_r^2 \\
        Corr\lbrace r(\vect x), r(\vect x')\rbrace \\
        &= exp\lbrace -\frac{\tau}{\xi_r}\rbrace
    \end{split}
\end{equation}

with $\tau = |\vect x - \vect x'|$

## Problem 3ac

Discretize the random field with  $L:\left[(1,30) \times (1,30)\right] \in D$ with parameters $\sigma_r^2$ and $\xi_r = 3$. And simulate realizations. 


\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig3a1-1} 

}

\caption{\label{fig:fig3a1} Realizations of field with $\xi = 3$, $\sigma^2 = 2$}\label{fig:fig3a1}
\end{figure}


Create four realizations that are presented in Figure \ref{fig:fig3a1}

## Problem 3bc
The theoretical variogram function can be written as: 
\begin{equation}
    \gamma_r(\tau) =  \sigma_r^2(1-\rho_r(\tau))
\end{equation}

We plot the emprical variogram from our realisations and compare to the theoretical variogram. 




\begin{figure}

{\centering \includegraphics{Exercise_1_files/figure-latex/fig3b1-1} 

}

\caption{\label{fig:fig3b1} Empirical variogram vs. theoretical variogram of realisation in Figure \ref{fig:fig3a1}.$\xi = 3$, $\sigma^2 = 2$}\label{fig:fig3b1}
\end{figure}

The results are plotted in Figure \ref{fig:fig3b1}. None of the empirical variograms seem to be very accurate. The seem to capture the shape of the theoretical variogram at low distances, however at large distances the empirical variograms vary alot, this is probably due to a low amount of points that have a large distance in between. 

\newpage 

## Problem 3de

Generate 36 point uniformly on $L$. And simulate a realisation over those 36 points. Will later take on these as exact observations. 




Use **likfit** to estimate the parameters of our field using our 36 exact observations.



The summary of the **likfit** is displayed inbelow:

```r
summary(likfit.36)
```

```
## Summary of the parameter estimation
## -----------------------------------
## Estimation method: maximum likelihood 
## 
## Parameters of the mean component (trend):
##   beta 
## 0.8442 
## 
## Parameters of the spatial component:
##    correlation function: exponential
##       (estimated) variance parameter sigmasq (partial sill) =  2.255
##       (estimated) cor. fct. parameter phi (range parameter)  =  1.869
##    anisotropy parameters:
##       (fixed) anisotropy angle = 0  ( 0 degrees )
##       (fixed) anisotropy ratio = 1
## 
## Parameter of the error component:
##       (estimated) nugget =  0
## 
## Transformation parameter:
##       (fixed) Box-Cox parameter = 1 (no transformation)
## 
## Practical Range with cor=0.05 for asymptotic range: 5.599689
## 
## Maximised Likelihood:
##    log.L n.params      AIC      BIC 
## "-64.19"      "4"  "136.4"  "142.7" 
## 
## non spatial model:
##    log.L n.params      AIC      BIC 
## "-67.12"      "2"  "138.2"  "141.4" 
## 
## Call:
## likfit(geodata = geo, ini.cov.pars = c(1, 1), cov.model = "exponential")
```

Using the 36 point we estimate $\hat\sigma^2 =$ 2.254746 and  $\hat \xi =$ 1.8692279.
We repeat the process in 3d for $n=6$, $n=64$, $n=100$. 


\newpage
The ML estimates of the parameters were: 

| n      | $\hat \sigma^2_{ML}$  | $\hat\xi_{ML}$      |
|--------|-----------------------|-------------------  |
| 9      | 1.1269907  |0.28203     |
| 36     | 2.254746 | 1.8692279   |
| 64     | 2.4823404 | 4.7979249   |
| 100    | 2.1700308|  3.3909212 |
|--------|-----------------------|-------------------  | 
| Actual | 2            | 3              | 
 

We see that the estimates becomes quite accurate at large values of $n$, and be very off at lower values of $n$

We plot the empirical variograms vs theoretical. 

![\label{fig:fig3d2} Variograms from uniform drawn points $\xi = 3$, $\sigma^2 = 2$, maximum likelihood fit](Exercise_1_files/figure-latex/unnamed-chunk-8-1.pdf) 

These variograms are displayed in Figure \ref{fig:fig3d2}. They seem to vary qutite much, but that can probably be blamed on a low amount of observation points. As we draw a very limited amount of points, some distances will be less represented than others, which will cause large variations in estimation. This can especially be seen for $n=9$


# Problem 3f
In this problem we have used maxmium likelihood to estiamte the parameters of the correlation funciton. We have found that using maximum likelihood to estiamte parameters seems to be relatively fast using the number of parameters and observations we have here (although **likfit** warns that the opposite might be the case). With enough data points maximum likelihood seems to be quite good at finding the correct parameters.



