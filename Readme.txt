
Our matlab code for paper "Low-rank Characteristic Tensor Density EstimationPart I: Foundations"


Effective non-parametric modeling and estimation of high-dimensional joint distributions is a fundamental problem in statistical learning and this paper reports basic results which can be applied in a  wide  range  of  applications mentioned in the introduction. The proposed probabilistic framework allows easy model marginalization, training using incomplete data, prediction of any variable from (a subset of) the other variables, and sampling, which can potentially benefit many research problems in applied statistics, data mining, and machine learning where density estimation can be used as a building block.
 
 
Our work shows that any smooth compactly supported multivariate PDF can be approximated by a finite tensor model. We do not assume a latent variable mixture model; instead, the latent variable factorization falls off from compactness of support and smoothness. This implication is a key result of our work, which was spelled out in Proposition 1 in our paper. 
	
 Our work further shows that under low-rank conditions, it is possible to infer the joint density of N>3 random variables from the joint density of triples of random variables (and, more generally, lower-dimensional distributions). We had shown this earlier  for discrete random variables, but in this paper we show that this remarkable Kolmogorov extension-type of result also holds for continuous random variables under certain reasonable conditions. 
	
We are the first to use tensor models for high-dimensional densities, where N is well above 3-10. Nobody believed that it is possible to use or compute low-rank tensor factorization in this regime. In fact, prior to our work, we do not know anyone who ever tried to decompose $100-$way tensors. Part of our novelty is that these models work, with remarkably low ranks, in these high dimensions (up to 256 in the paper).

