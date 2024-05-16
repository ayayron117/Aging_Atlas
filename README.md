# *c.elegan* Aging Atlas

# [Bootstrapping](https://github.com/ayayron117/Aging_Atlas/blob/main/bootstrap.md)

<div align='justify'>
For each tissue of each genotype of each time point, 15 cells were randomly sampled and the counts associated with each of the 26208 genes were aggregated. This was done 101 times and each result was stored as a column in a matrix, where each row was a gene. I used the sample() function which used the Mersenne-Twister method for random number generation, the default method used by version of R that I was using at the time (4.2.1). The cells were sampled without replacement.
</div>


