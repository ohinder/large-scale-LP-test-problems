Repository for generating large-scale LP test instances for PDLP MPC paper.

More examples and documentation will be added in the near future.

# Creating the test instances

To create small versions of the problem instances (for testing) run:

```shell
$ ./create_small_instances.sh
```

To create large versions of the problem instances (for testing) run. This is going to take a while and you will need a machine 
with >= 32GB of RAM.

```shell
$ ./create_large_instances.sh
```

# Background on problems

## Synthetic model of the supply chain for a large retailer

### Motivation

This problem is motivated by optimizing the supply chain of a large retailer
over relatively moderate time period (e.g., a month).
The retailer sells $K$ different products (commodities) which 
it purchases from many different factories.
The retailer has $F$ factories for each commodity $k = 1, \dots, K$
and $W$ total warehouses.
The retailer needs to meet demand at its $S$ different stores
shipping the goods from the factories to the warehouses
and then to the stores.
It is assumed each commodity has a known demand at each store.
There is a shipping cost for transporting the goods from the factories to
the warehouses, and from the warehouses to the stores which is 
proportional to the distance the goods are transported.
Each commodity is produced at $F$ different factories.
The locations of these factories are different for each
commodity. It is imagined that the retailer 
is operating during peak season where each of their warehouses
has a limited normal processing capacity of goods (the average of
which is less than total demand).
However, the number of goods processed by a warehouse can be increased
by the use overtime (but the retailer must incurs extracost by doing so).
The retailer wishes to meet demand while minimizing total costs.
We allow fractional flows of goods.
This is a (fractional) multicommodity flow problem. 

We have tried to make this as simple as possible while maintaining
some realisism. In practice, one can imagine more complicated 
and larger models incorporating, for example,
the multiperiod aspect of the problem or the issue of deciding where 
to build new warehouses.

### Decision variables

Let $u_{k,f,w} \ge 0$ be the flow from factory $(k,f)$ to warehouse $w$ for commodity $k$.
Let $v_{k,w,s} \ge 0$ be the flow from warehouse $w$ to store $s$ for commodity $k$.
Let $x_{w} \ge 0$ be the amount of overtime at warehouse $w$.

### Parameters

Let $g_{k,f}$ be the location of factory $f$ for commodity $k$.
Let $h_{w}$ be the location of warehouse $w$.
Let $q_{s}$ be the location of store $s$.
The values are generated by sampling randomly
from a uniform distribution between zero and one.

Let $d_{k,s}$ be the demand for store $s$ for commodity $k$.
The values for $d_{k,s}$ are generated in a two step process.
If the first step we generate average demands for each 
commodity
$a_{k} = 100 \exp(Z_k)$ where $Z_1, \dots, Z_K$ are sampled indepedently from a standard normal distributation.
Then demands at each store $s$ for commodity $k$ are realized by setting
$d_{k,s} \sim \text{Possion}(a_k)$.
This two step process generates more realistic data as some products have very high demand while most products have very low demand.

Let $m_{k,f}$ be the maximum production of commodity $k$ at factory
$f$. This is generated in the code using 
the formula $m_{k,f} = \frac{1}{F} \sum_{k=1}^K \sum_{s=1}^S d_{k,s}$ such 
that the factories produce the same quantity of each
code and total demand equals total supply.

Let $\gamma$ be the normal processing capacity of a
warehouse. This is generated in the code using 
the formula $\gamma = \frac{0.95}{W} \times \sum_{k=1}^K \sum_{s=1}^S d_{k,s}$
so that without using any overtime the warehouses could meet $95\%$
of total demand.
Let $\theta$ be the cost of additional overtime.
By default this is set to be $0.3$.

### Optimization model

Minimize shipping costs plus overtime costs:
$$
\min \sum_{k=1}^K \sum_{f=1}^F \sum_{w=1}^W \| g_{k,f} - h_{w} \|_2 u_{k,f,w} + \sum_{s=1}^S \sum_{w=1}^W \| h_{w} - q_s \|_2 u_{k,f,w} + \theta \sum_{w=1}^W x_{w}.
$$
Flow from each factory $(f,k)$ does not exceed supply:
$$
\sum_{w=1}^W u_{k,f,w} \le m_{k, f}.
$$
For each warehouse $w$, if the normal operating capacity $\gamma$ is exceeded then we
must use overtime:
$$
\sum_{k=1}^K \sum_{f=1}^F u_{k,f,w} \le \gamma + x_{w}.
$$
For each commodity $k$, flow into warehouses $w$ equals flow out of warehouse $w$:
$$
\sum_{f=1}^F u_{k,f,w} = \sum_{s=1}^S v_{k,w,s}.
$$
Demand for each commodity $k$ at store $s$ is met:
$$
\sum_{w=1}^W v_{k,w,s} \ge d_{k,s}.
$$

### Generating plots of the optimal solution

For small instances you can also solve the problem using HiGHS and generate
plots of the optimal solution for a random selection of the commodities:

```shell
$ julia-1.7 generate-multicommodity-flow.jl \
    --output_file problem-instances/multicommodity-flow-tiny-test-instance.mps.gz \
    --optimize_model \
    --folder_for_plots plot_optimal_flows
```

This was used during the development process to verify the model was producing reasonable optimal solutions.

## Locating heat sources from a small number measurements

### Motivation

Consider a homogenous cube $[0,1]^3$ of material with constant temperature on its boundary.
Within the material several heat sources have been randomly located.
Using a small number of temperature measurements in the material 
from predefined locations we wish to recover the location
of these heat sources. To make the problem the problem more tractable
there is only a small number of possible locations for these heat sources.
We assume the material is at equilibrium.
If locations of the heat sources are known (and the size of their input)
then we can solve Possion's equation to calculate the temperature profile
through the material:
$$
\frac{d^2 u}{d x^2} + \frac{d^2 u}{d y^2} + \frac{d^2 u}{d z^2} = -q(x, y, z)
$$
where $q$ represents the heat source inputs and $u$ is the temperature profile.
This is an inverse problem: we want to recover the heat sources from the 
observations. To solve this problem we will discretize the PDE into
$N \times N \times N$ grid points.

### Decision variables

The temperature of the material at location $(i / N,j / N, k / N)$ is given by
$u_{i,j,k}$. The heat injected at location $(i / N,j / N, k / N)$ is
given by $q_{i,j,k} \ge 0$.

### Optimization model

As there is relatively few heat sources we minimize the sum of 
$q_{i,j,k}$, this naturally makes the model sparse:
$$
\min \sum_{i,j,k} q_{i,j,k}
$$
We then build a finite element approximation of Possion's equation:
$$
\frac{u_{i+1,j,k} - 2 u_{i,j,k} + u_{i-1,j,k}}{h^2} + \frac{u_{i,j+1,k} - 2 u_{i,j,k} + u_{i,j-1,k}}{h^2} + \frac{u_{i,j,k+1} - 2 u_{i,j,k} + u_{i,j,k-1}}{h^2} = -q_{i,j,k}
$$
and setup the boundary conditions
$$
u_{1,j,k} = 0 \quad u_{i,1,k} = 0 \quad u_{i,j,1} = 0
$$
$$
u_{N,j,k} = 0 \quad u_{i,N,k} = 0 \quad u_{i,j,N} = 0.
$$
Then at the measurement locations $M$ we have
$$
u_{i,j,k}^\star = u_{i,j,k}  \quad \forall (i,j,k) \in M
$$

### Instance generation

We generate $v_1, \dots, v_m$ measurement locations uniformly at random
from $[0,1]^3$ and potential locations for the heat sources $w_1, \dots, w_k$ uniformly at random
from $[0,1]^3$. Then
from the $k$ potential locations for the heat sources we choose uniformly at random $l$ location different to place the heat sources, denoted by the set $\mathcal{L}$ where $l < k$. At each location with a heat source we generate the 
total inputted heat uniformly between zero and one.
We then generate our grid and round the positional vectors to the nearest point in the grid. 
This yields a set $M$ of indicies for the measurement locations.
Finally, we solve the discretized Possion's equation to calculate the true temperature distribution $u^\star$.

## Statistical matching with covariate balancing constraints for making causal inference on observational data

### Motivation

In observational studies we are often wish to perform causal inference where we try to understand the 
impact of particular treatment on patient outcomes. However, one difficulty is that the patients
that recieve the treatment can be a very different group of patients from those that did not recieve
the treatement. For example, perhaps there is a confounding variable -- how sick the patient is
which affects whether they recieve the treatement and their outcome.  
Matching is a popular approach in statistics where treated patients are matched to control patients
with similar covariate values. This produces a new smaller sample where the patients in both the
treatment and control are similar -- resembling a randomized experiment. 
However, one problem with this approach is that the covariates in the treatment group may not match covariates in the subsampled control group [C].
For this reason, it has been proposed that ones finds the optimal matching subject
to additional constraints that the in the covariates in the subsampled control group match the treatment [A,B]. 
This can be modelled as integer program [A] but in practice the linear programming
relaxation can be used to quickly generate matchings for large
datasets [B].
We generate synthetic versions of this problem for testing.

### Parameters

There are $i = 1, \dots, n$ samples in the treatment group and $j = 1, \dots, m$ samples in the control group. We assume that $m > n$. Moroever, there are $k=1,\dots,d$ covariate values for each sample.

$A_{ik}$ for $i=1,\dots,n$ and $k=1,\dots,d$ is the $k\text{th}$ covariate value for sample $i$ of the treatment group.

$B_{ik}$ for $i=1,\dots,n$ and $k=1,\dots,d$ is the $k\text{th}$ covariate value for sample $j$ of the control group.

$q$ is the total number of samples to be included in the subsample.

Furthermore, define the first moment of the original sample as
$$
\bar{A}_k = \frac{1}{n+m} \left(\sum_{i=1}^n A_{ik} + \sum_{j=1}^m B_{jk} \right)
$$
Also, define the second moment of the orignal sample as
$$
M_{kl} = \frac{1}{n+m} \left( \sum_{i=1}^n A_{ik} A_{il} + \sum_{j=1}^m B_{ik} B_{il} \right)
$$
Let $E$ be the set of edges, i.e., the pairs of treatment and control
samples that we allow to be matched.

### Decision variables

$x_{ij}$ is one if sample $i$ from the treatment group is matched to sample $j$ from the control group and zero otherwise.

$w_{j}$ for $j = 1, \dots, n$ is one if sample $j$ from the control group is matched and zero otherwise.

To turn it into a linear program, we relax the integrality requirements and only require that $x_{ij}$ and $w_j$ are in $[0,1]$.
This is a popular approach for solving large-scale versions of these problems [B].

### Optimization model

Minimize the total weight of the assignment:
$$
\sum_{(i,j) \in E} \| A_{i\cdot} - B_{i\cdot} \| x_{ij}
$$
Matching constraints for treatment group:
$$
\sum_{j=1}^{m} x_{ij} = 1 \quad i = 1, \dots, n
$$
Matching constraints for control group:
$$
\sum_{i=1}^{n} x_{ij} = w_{j} \quad i = 1, \dots, m
$$
Covariate first moments of subsampled control approximately match the treatment:
$$
-\epsilon \le  \frac{1}{n} \sum_{j=1}^m B_{jk} w_j - \bar{A}_k \le \epsilon \quad k = 1, \dots, d
$$
Covariate second moments approximately match original sample:
$$
-\epsilon \le  \frac{1}{n} \sum_{j=1}^m B_{jk} B_{jl} w_j - M_{kl} \le \epsilon \quad k = 1, \dots, d \quad l = 1, \dots, d
$$

### Instance generation

We generate $A_{ik}$ from a standard normal distribution.
We then create a shift vector:
$$
v_{k} \sim N(0,0.1) \quad k = 1, \dots, d
$$
and then set
$$
B_{ik} \sim N(v_k, 1).
$$
This shift vector ensures that the covariate values of $A$ and $B$ have different distributions.
Finally, to choose the edges we use a k-d tree to find 
the $t$ closest control samples to each treatment sample (where $t$ is say $10$).
This reduces the number of edges we include in the problem, 
allowing us to focus on the most promising matches.
This make the problem size dramatically smaller but, 
in our experience,
has minimal impact on the optimal objective.

### References

[A] Zubizarreta, José R. "Using mixed integer programming for matching in an observational study of kidney failure after surgery." Journal of the American Statistical Association 107.500 (2012): 1360-1371.

[B] https://cran.r-project.org/web/packages/designmatch/index.html

[C] Ali, M.S., Groenwold, R.H., Belitser, S.V., Pestman, W.R., Hoes, A.W., Roes, K.C., de Boer, A. and Klungel, O.H., 2015. Reporting of covariate selection and balance assessment in propensity score analysis is suboptimal: a systematic review. Journal of clinical epidemiology, 68(2), pp.122-131.