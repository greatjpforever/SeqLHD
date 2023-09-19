# Sequentially refined Latin hypercube designs with flexibly and adaptively chosen sample sizes

This project includes 'main_compare_rmse.m' for conducting numerical experiments to assess the quality of various sequential designs, as well as multiple algorithms for generating Latin Hypercube Designs (LHD).
These algorithms can help you generate uniformly distributed sample points, which are useful for experimental design, parameter optimization, and computer simulations.
Let's now shift our focus to these specific algorithms.

## Algorithm List

1. `seqlhd_alg1`: Initial Latin Hypercube Design algorithm.

2. `seqlhd_alg2`: Extended Latin Hypercube Design algorithm.

3. `seqlhd_alg3`: Another extended Latin Hypercube Design algorithm.

## Parameter

- Input arguments of seqlhd_alg1:
  - n1: Positive integer indicating the number of runs.
  - d: Positive integer indicating the number of dimensions.
  - m1: (Optional) Positive integer indicating the distance parameter. Default: ceil(n1/5)+2

- Input arguments of seqlhd_alg2:
  - E1: n1-by-d matrix representing the original Latin hypercube design.
  - n2: Positive integer indicating the number of additional samples to be added
  - m1: Positive integer indicating the distance parameter of E1. Default: ceil(n1/5)+2

- Input arguments of seqlhd_alg3:
  - E2: N2-by-d matrix representing the original Latin hypercube design.
  - n3: Positive integer indicating the number of additional samples to be added.

Note: Pay attention that there are some additional requirements regarding the parameters depending on how you use these algorithms. These requirements are detailed in the usage and examples provided below.

## Usage and example

### generate an initial LHD

```matlab
% Generate a 2-dimensional Latin Hypercube Design matrix with 3 sample points per dimension
% first way, use the default m1
n1 = 3;
d = 2;
D1 = seqlhd_alg1(n1, d); % m1 will default to ceil(n1/5)+2
disp(D1)
>>
    0.6562    0.8967
    0.8885    0.3080
    0.0276    0.6342
% second way, use other m1
n1 = 3;
d = 2;
m1 = 4;
D1 = seqlhd_alg1(n1, d, m1); 
disp(D1)
>> 
    0.3628    0.6847
    0.9055    0.3982
    0.0619    0.1280
```

### Sequentially generate a LHD with three stages

In this way, we use seqlhd_alg1, seqlhd_alg2 and seqlhd_alg3 to generate a three-stage Latin hypercube design;
$n_2$ needs to be greater than $m_1$.
Moreover, if $1/n_1+1/(n_1+m_1)\geq 3/(n_1+n_2)$, then $n_3$ needs to be greater than $n_1+n_2$.
Otherwise, let $d=\min\{1/(2N_2), n_2(1-h_1 N_2)/(n_1N_2)\}$. Then $n_3$ needs to be greater than $1/d-n_1-n_2$.

```matlab
% generate the first stage using seqlhd_alg1
n1 = 2;
d = 2;
m1 = ceil(n1/5)+2;
E1 = seqlhd_alg1(n1, d, m1);
disp(E1);
>>
    0.4546    0.2800
    0.9102    0.6921

% extend E1 using seqlhd_alg2
n2 = 3; % n2 need be greater than m1 
E2 = seqlhd_alg2(E1, n2, m1);
disp(E2);
>>
    0.4546    0.2800
    0.9102    0.6921
    0.3358    0.0405
    0.7641    0.4324
    0.1353    0.9372

% extend E2 using seqlhd_alg3
% n3 need be greater than n1+n2, because 1/n1+1/(n1+m2)>3/(n1+n2)
n3 = 19; 
E3 = seqlhd_alg3(E2,n3)
disp(E3);
    0.4546    0.2800
    0.9102    0.6921
    0.3358    0.0405
    0.7641    0.4324
    0.1353    0.9372
    0.0193    0.4074
    0.2361    0.2425
    0.9510    0.1222
    0.2615    0.6101
    0.5008    0.9009
    0.6763    0.5720
    0.6238    0.7921
    0.4546    0.2800
    0.9102    0.6921
    0.3358    0.0405
    0.7641    0.4324
    0.1353    0.9372
    0.0193    0.4074
    0.2361    0.2425
    0.9510    0.1222
    0.2615    0.6101
    0.5008    0.9009
    0.6763    0.5720
    0.6238    0.7921

```

### Alternative use: Sequentially generate a LHD with two stages

**In this way, we use seqlhd_alg1 and seqlhd_alg3 to generate a two-stage Latin hypercube design;**
**n2 needs to be greater than m1;**

```matlab
% generate the first stage using seqlhd_alg1
n1 = 3;
d = 2;
E1 = seqlhd_alg1(n1, d); % m1=ceil(n1/5)+2=2
disp(E1);
>>
    0.7041    0.9146
    0.4799    0.1135
    0.1595    0.3777

% generate the second stage using seqlhd_alg3 
n2 = 4; % n2 need be greater than m1
E2 = seqlhd_alg3(E1, n2);
disp(E2);
>>
    0.7041    0.9146
    0.4799    0.1135
    0.1595    0.3777
    0.0811    0.6390
    0.3107    0.7976
    0.9810    0.5144
    0.8452    0.2771
```

### Notes

- Please adjust the input parameters and example data according to your needs.
- Input parameters should conform to the requirements of the functions, such as positive integers or non-negative integers.
- Make sure that the MATLAB or Octave environment is properly set up before calling the algorithm functions.
- Ensure that the algorithm code is correctly copied into your project directory and properly referenced in your code.
- Adjust the input parameters and example data according to your specific requirements to generate a Latin Hypercube Design matrix suitable for your project.
- Contact the Author
For any questions or suggestions, please contact the author:
Email: <gongjunpeng21@mails.ucas.ac.cn>
- Feel free to modify and enhance the README example to better suit your project's needs.

Thank you for using this project! If you find it helpful, please provide your support and feedback.
