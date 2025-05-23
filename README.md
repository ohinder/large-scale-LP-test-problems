Repository for generating large-scale linear programming (LP) test instances for the PDLP MPC paper.

# One-time setup
The code was tested using Julia 1.10.0. To instantiate all the required packages execute 
``` shell
$ julia --project=. -e "import Pkg; Pkg.instantiate()"
```

# Creating the test instances

To create small versions of the problem instances (for testing) run:

```shell
$ ./create_small_instances.sh
```

To create large versions of the problem instances run:

```shell
$ ./create_large_instances.sh
```

This is going to take a while and you will need a machine 
with $`\ge 512`$GB of RAM.

You can modify the instance generation properties and seeds.
Take a look at the above shell files for examples of parameters
that you can pass to the generators. For more information
on the parameters for each generator you can use `--help`,
for example, for `generate-multicommodity-flow.jl` you can run

```shell
julia --project generate-multicommodity-flow.jl --help
```

# Background on problems

## Synthetic model of the supply chain for a large retailer

### Motivation

This problem is motivated by optimizing the supply chain of a large retailer
over a relatively moderate time period (e.g., a month). The retailer sells $K$
different products (commodities), which it purchases from many different
factories. The retailer has $F$ factories for each commodity $k = 1, \dots, K$,
and $W$ total warehouses. The retailer needs to meet demand at its $S$ different
stores, shipping the goods from the factories to the warehouses and then to the
stores. It is assumed each commodity has a known demand at each store. There is
a shipping cost for transporting the goods from the factories to the warehouses,
and from the warehouses to the stores, which is proportional to the distance the
goods are transported. Each commodity is produced at $F$ different factories.
The locations of these factories are different for each commodity. It is assumed
that the retailer is operating during peak season where each of their warehouses
has a limited normal processing capacity of goods (the average of which is less
than total demand). However, the number of goods processed by a warehouse can be
increased by the use over time (but the retailer incurs in an extra cost by
doing so). The retailer wishes to meet demand while minimizing total costs. We
allow fractional flows of goods, so this is a (fractional) multicommodity flow
problem.

We have tried to make this as simple as possible while maintaining
some realism. In practice, one can imagine more complicated 
and larger models incorporating, for example,
the multiperiod aspect of the problem or the issue of deciding where 
to build new warehouses.

For the `create_large_instances.sh` we use 20000 commodities, 100 warehouses and 1000 stores. For context, Walmart has over 4000 stores in the United States, and the company has over 200 distribution centers [E]. Moreover, the average Walmart
store sells around 120,000 different products [D].

### Decision variables

Let $u_{k,f,w} \ge 0$ be the flow from factory $(k,f)$ to warehouse $w$ for commodity $k$.
Let $v_{k,w,s} \ge 0$ be the flow from warehouse $w$ to store $s$ for commodity $k$.
Let $x_{w} \ge 0$ be the amount of overtime at warehouse $w$.

### Parameters

The locations below belong to the 2D square $[0, 1] \times [0, 1].$
Let $g_{k,f}$ be the location of factory $f$ for commodity $k$.
Let $h_{w}$ be the location of warehouse $w$.
Let $q_{s}$ be the location of store $s$.
The values of these parameters are generated by sampling each component randomly from a uniform distribution between zero and one.
 
Let $d_{k,s}$ be the demand for store $s$ for commodity $k$.
The values for $d_{k,s}$ are generated in a two step process.
If the first step we generate average demands for each 
commodity
$a_{k} = 100 \exp(Z_k)$ where $Z_1, \dots, Z_K$ are sampled indepedently from a standard normal distribution.
Then, demands at each store $s$ for commodity $k$ are realized by setting
$d_{k,s} \sim \text{Possion}(a_k)$.
This two step process generates more realistic data as some products have very high demand while most products have very low demand.

Let $m_{k,f}$ be the maximum production of commodity $k$ at factory
$f$. This is generated in the code using 
the formula $`m_{k,f} = \frac{1}{F} \sum\limits_{k=1}^K \sum\limits_{s=1}^S d_{k,s}`$ such 
that the factories produce the same quantity of each
product and total demand equals total supply.

Let $\gamma$ be the normal processing capacity of a
warehouse. This is generated in the code using 
the formula $`\gamma = \frac{0.95}{W} \times \sum\limits_{k=1}^K \sum\limits_{s=1}^S d_{k,s}`$
so that without using any overtime the warehouses could meet $95\%$
of total demand.
Let $\theta$ be the cost of additional overtime.
By default this is set to be $0.3$.

### Optimization model

Minimize shipping costs plus overtime costs:

```math
\min \sum_{k=1}^K \sum_{f=1}^F \sum_{w=1}^W \| g_{k,f} - h_{w} \|_2 u_{k,f,w} + \sum_{s=1}^S \sum_{w=1}^W \| h_{w} - q_s \|_2 u_{k,f,w} + \theta \sum_{w=1}^W x_{w}.
```
With the following constratints. Flow from each factory $(f,k)$ does not exceed supply:
``` math
\sum_{w=1}^W u_{k,f,w} \le m_{k, f}.
```
For each warehouse $w$, if the normal operating capacity $\gamma$ is exceeded then we
must use overtime:
``` math
\sum_{k=1}^K \sum_{f=1}^F u_{k,f,w} \le \gamma + x_{w}.
```
For each commodity $k$, flow into warehouses $w$ equals flow out of warehouse $w$:
``` math
\sum_{f=1}^F u_{k,f,w} = \sum_{s=1}^S v_{k,w,s}.
```
Demand for each commodity $k$ at store $s$ is met:
``` math
\sum_{w=1}^W v_{k,w,s} \ge d_{k,s}.
```

### Generating plots of the optimal solution

For small instances you can also solve the problem using HiGHS and generate
plots of the optimal solution for a random selection of the commodities:

```shell
$ julia --project generate-multicommodity-flow.jl \
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
of these heat sources. To make the problem more tractable
there is only a small number of possible locations for these heat sources.
We assume the material is at equilibrium.
If locations of the heat sources are known (and the size of their input)
then we can solve Possion's equation to calculate the temperature profile
through the material:
``` math
\frac{d^2 u}{d x^2} + \frac{d^2 u}{d y^2} + \frac{d^2 u}{d z^2} = -q(x, y, z)
```
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
``` math
\min \sum_{i,j,k} q_{i,j,k}
```
We then build a finite element approximation of Possion's equation:
``` math
\frac{u_{i+1,j,k} - 2 u_{i,j,k} + u_{i-1,j,k}}{h^2} + \frac{u_{i,j+1,k} - 2 u_{i,j,k} + u_{i,j-1,k}}{h^2} + \frac{u_{i,j,k+1} - 2 u_{i,j,k} + u_{i,j,k-1}}{h^2} = -q_{i,j,k}
```
and setup the boundary conditions

``` math
u_{1,j,k} = 0 \quad u_{i,1,k} = 0 \quad u_{i,j,1} = 0
```
``` math
u_{N,j,k} = 0 \quad u_{i,N,k} = 0 \quad u_{i,j,N} = 0.
```
Then at the measurement locations $M$ we have
``` math
u_{i,j,k}^\star = u_{i,j,k}  \quad ~~\quad\forall (i,j,k) \in M.
```

### Instance generation

We generate $v_1, \dots, v_m$ measurement locations uniformly at random
from $[0,1]^3$ and potential locations for the heat sources $w_1, \dots, w_n$ uniformly at random
from $[0,1]^3$. Then
from the $n$ potential locations for the heat sources we choose uniformly at random $l$ location different to place the heat sources, denoted by the set ${L}$ where $l < n$. At each location with a heat source we generate the 
total inputted heat uniformly between zero and one.
We then generate our grid and round the positional vectors to the nearest point in the grid. 
This yields a set $M$ of indices for the measurement locations.
Finally, we solve the discretized Possion's equation to calculate the true temperature distribution $u^\star$.

### Validation of recovery

One can also validate recovery of the true solution by using the script `validate-heat-source-solution.jl`.
This script compares the temperature profile $u$ of the true solution (computed during instance generation)
with the solution obtained by PDLP from solving the optimization problem.

## Statistical matching with covariate balancing constraints for making causal inference on observational data

### Motivation

In observational studies we often wish to perform causal inference where we try to understand the 
impact of particular treatment on patient outcomes. However, one difficulty is that the patients
that receive the treatment can be a very different group of patients from those that did not receive
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

$B_{ij}$ for $j=1,\dots,m$ and $k=1,\dots,d$ is the $k\text{th}$ covariate value for sample $j$ of the control group.

$q$ is the total number of samples to be included in the subsample.

Furthermore, define the first moment of the treatment sample as
``` math
\bar{A}_k = \frac{1}{n} \left(\sum_{i=1}^n A_{ik} \right)
```
Also, define the second moment of the treatment sample as
``` math
M_{kl} = \frac{1}{n} \left( \sum_{i=1}^n A_{ik} A_{il} \right)
```
We consider a bipartite graph with vertices given by $V = V_t \cup V_c$, where $V_t \cap V_c \cap =\emptyset$, and $V_t$ and $V_c$ denote paients in the treatment and control groups, respectively. We define a subset of the edges that can be used for matching $E \subset \{\{u, v\} \mid u \in V_t \text{ and } v \in V_c\}$; see [Instance generation](### Instance generation) below for more details about how this subset is selected.


### Decision variables

$x_{ij}$ is one if sample $i$ from the treatment group is matched to sample $j$ from the control group and zero otherwise.

$w_{j}$ for $j = 1, \dots, n$ is one if sample $j$ from the control group is matched and zero otherwise.

To turn it into a linear program, we relax the integrality requirements and only require that $x_{ij}$ and $w_j$ are in $[0,1]$.
This is a popular approach for solving large-scale versions of these problems [B].

### Optimization model

Minimize the total weight of the assignment:
``` math
\sum_{(i,j) \in E} \| A_{i\cdot} - B_{i\cdot} \| x_{ij}.
```
Subject to the following constraints. Matching constraints for treatment group:
``` math
\sum_{j=1}^{m} x_{ij} = 1 \quad i = 1, \dots, n.
```
Matching constraints for control group:
``` math
\sum_{i=1}^{n} x_{ij} = w_{j} \quad i = 1, \dots, m.
```
Covariate first moments of subsampled control approximately match the treatment:
``` math
-\epsilon \le  \frac{1}{n} \sum_{j=1}^m B_{jk} w_j - \bar{A}_k \le \epsilon \quad k = 1, \dots, d.
```
Covariate second moments approximately match original sample:
``` math
-\epsilon \le  \frac{1}{n} \sum_{j=1}^m B_{jk} B_{jl} w_j - M_{kl} \le \epsilon \quad k = 1, \dots, d \quad l = 1, \dots, d.
```

### Instance generation
We generate $A_{ik}$ from a standard normal distribution.
We then create a shift vector:
``` math
v_{k} \sim N(0,0.1) \quad k = 1, \dots, d
```
and then set
``` math
B_{ik} \sim N(v_k, 1).
``` 
This shift vector ensures that the covariate values of $A$ and $B$ have different distributions.
Finally, to choose the edges we use a [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) to find 
the $t$ closest control samples to each treatment sample (where $t$ is set to $10$).
This reduces the number of edges we include in the problem, 
allowing us to focus on the most promising matches.
This make the problem size dramatically smaller but, 
in our experience,
has minimal impact on the optimal objective.

## Robust Production Inventory Problem

### Motivation
We consider the classic robust production-inventory problems from [F]. The class of problems considers a firm with a central warehouse and $E$ factories that aim to satisfy uncertain demand for a single product over a selling season.  The selling season of the firm's product is discretized into $T$ time periods, which are spaced equally over the selling season. 

In each time period,  the firm sequentially performs the following three steps: 

1. The firm replenishes the inventory level at the central warehouse by producing additional products at their  factories.  Each  factory produces the additional units with  zero  lead time, and the additional units  are stored  immediately in the central warehouse. 

2. The firm observes the customer demand  at the central warehouse. 

3. The firm verifies that the remaining  inventory  in the warehouse lies within a  pre-specified interval.

In addition to satisfying the constraints on the inventory level in the central warehouse at the end of each time period, the firm's production decisions must be in a certain range. 

The goal of the firm is to satisfy the customer demand at minimal cost while satisfying production and warehouse constraints. 

### Parameters

Let  $x_{te} \ge 0$ denote the number of product units that the firm decides to produce at each of the factories $e = 1, ..., E$ at a per-unit cost of $c_{te}$. 

The demand at the central warehouse is denoted by 
``` math
\zeta_{t+1} \in {U}_{t+1} \equiv [\underline{D}_{t+1},\bar{D}_{t+1}]
```
which must be satisfied immediately without backlogging from the inventory in  the central warehouse. The lower and upper bounds in the uncertainty set, denoted by  $\underline{D}\_{t+1} < \bar{D}\_{t+1},$ capture the minimum and maximum level of customer demand that the firm anticipates  receiving in each time period $t$.

The remaining inventory level in the central warehouse at the end of each time period $t = 1, ..., T$ must satisfy
``` math
    V_{\text{min}} \le  v_1 + \sum_{\ell =1}^t  \sum_{e=1}^E x_{\ell e} - \sum_{s=2}^{t+1} \zeta_s \le V_{\text{max}},
```
where $v_1$ is the initial inventory level in the central warehouse at the beginning of the selling horizon, the second term is the cumulative number of product units that have been produced at the factories up through time period $t$, and the third term is the cumulative customer demand that has been observed at the central warehouse up through time period $t$.

### Decision variables
We utilize a linear decision rule for the above robust optimization problems. More specifically, we set 
``` math
x_{te}=\sum_{s=1}^t y_{t,s,e} \zeta_s.
```
 The decision variable is the linear decision rule parameter $y_{t,s,e}$.

### Optimization model
Minimize the worst-case cost with uncertainty:
``` math
\underset{\substack{y_{t,1},\ldots,y_{t,t} \in R^E:\;  ~~\quad\forall t = 1, ..., T}}{\textnormal{minimize}} \max_{\zeta_1 \in {U}_1,\ldots,\zeta_{T+1} \in {U}_{T+1}} \left \{ \sum_{t=1}^T \sum_{e=1}^E c_{te} \left( \sum_{s=1}^t    y_{t,s,e} \zeta_s\right) \right \}.
```
Subject to the following constraints. Maximal total production level for each factory:
``` math
\sum_{t=1}^T  \left( \sum_{s=1}^t y_{t,s,e} \zeta_s \right) \le  Q_e   ~~\quad\forall e = 1, ..., E.
```
Maximal and minimal production level for each factor at a time period:
``` math
0 \le \left( \sum_{s=1}^t y_{t,s,e} \zeta_s \right) \le  p_{te}  ~~\quad\forall  e = 1, ..., E,\; t = 1, ..., T.
```
The remaining  inventory  in the warehouse lies within a  pre-specified interval:
``` math
V_{\textnormal{min}} \le  v_1 + \sum_{\ell=1}^t \sum_{e=1}^E \left( \sum_{s=1}^\ell y_{\ell,s,e} \zeta_s \right)  - \sum_{s=2}^{t+1} \zeta_s  \le  V_{\textnormal{max}} ~~\quad\forall  t = 1, ..., T.
```
### Instance generation
We generate the instance following [G], which generalized those from [F]. More specifically, we generate instances in which the customer demand and production costs of a new product follow a cyclic  pattern due to seasonality over a selling horizon of one year. 

Given a discretization of the selling season into $T$ stages, the  customer demand in 
``` math
    \phi_t = 1 + 0.5 \sin\left(\frac{2 \pi (t-2)}{T} \right), \theta_t = 0.2,  {U}_t =  \left[ \frac{1000  (1 - \theta)  \phi_t}{T / 24}, \frac{1000 (1+\theta)  \phi_t}{T/24}  \right], 
```

Given $E$ factories available to the firm, the production costs and capacities for each stage $t = 1, ..., T$ and each factory $e = 1, ..., E$ are 
``` math
    c_{te} = \left(1 + \frac{e-1}{E-1} \right) \phi_t,  p_{te} = \frac{567}{\left( T/24\right)\left(E/3 \right)},  Q_e = \frac{13600}{E / 3},
```
and the capacities and initial inventory at the central warehouse are 
``` math
V_{\text{min}} = 500, \quad V_{\text{max}} = 2000, \quad v_1 = 500.
```

### References

[A] Zubizarreta, José R. "Using mixed integer programming for matching in an observational study of kidney failure after surgery." Journal of the American Statistical Association 107.500 (2012): 1360-1371.

[B] https://cran.r-project.org/web/packages/designmatch/index.html

[C] Ali, M.S., Groenwold, R.H., Belitser, S.V., Pestman, W.R., Hoes, A.W., Roes, K.C., de Boer, A. and Klungel, O.H., 2015. Reporting of covariate selection and balance assessment in propensity score analysis is suboptimal: a systematic review. Journal of clinical epidemiology, 68(2), pp.122-131.

[D] https://www.zippia.com/advice/walmart-statistics/#:~:text=every%20single%20day.-,In%20fact%2C%20just%20one%20Walmart%20Supercenter%20in%20the%20U.S.%20serves,items%20available%20on%20the%20shelves.

[E] https://corporate.walmart.com/news/2022/06/03/a-new-era-of-fulfillment-introducing-walmarts-next-generation-fulfillment-centers

[F] Aharon Ben-Tal, Alexander Goryashko, Elana Guslitzer, and Arkadi Nemirovski. Adjustable robust
solutions of uncertain linear programs. Mathematical Programming, 99(2):351–376, 2004.

[G] Lu, Haihao, and Brad Sturt. "On the Sparsity of Optimal Linear Decision Rules in Robust Inventory Management." arXiv preprint arXiv:2203.10661 (2022).
