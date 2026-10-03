set.seed(1)
out_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"

###Set parameters for simulations
n <- 50000
p <- 10000
#prob <- runif(p, 0.1, 0.9)
prob <- rep(0.5, p)
a <- rep(0, p)
d <- rep(0, p)
h2 <- 0.5

A1A1 <- 2
A1A2 <- 1
A2A2 <- 0

###Simulate genotypes
SNPs <- matrix(NA, n, p)

for(j in 1:p){
  SNPs[, j] <- rbinom(n=n, size=2, prob=prob[j])
}

pheno <- rnorm(n, 0, 1)

###Plot genetic values and phenotypes
a_hat <- vector("numeric", p)
d_hat <- vector("numeric", p)

for(j in 1:p){
  a_hat[j] <- (mean(pheno[which(SNPs[, j] == A1A1)]) - mean(pheno[which(SNPs[, j] == A2A2)]))/2
  d_hat[j] <- mean(pheno[which(SNPs[, j] == A1A2)]) - (mean(pheno[which(SNPs[, j] == A1A1)]) + mean(pheno[which(SNPs[, j] == A2A2)]))/2
}  

counts <- apply(SNPs, 2, table)
n_A1A1 <- counts[3, ]
n_A1A2 <- counts[2, ] 
n_A2A2 <- counts[1, ]

var_A1A1 <- numeric(p)  # Variance for A1A1 genotype
var_A1A2 <- numeric(p)  # Variance for A1A2 genotype
var_A2A2 <- numeric(p)  # Variance for A2A2 genotype

for (j in 1:p) {
  # Subset phenotypes by genotype for locus j
  pheno_A1A1 <- pheno[SNPs[, j] == A1A1]
  pheno_A1A2 <- pheno[SNPs[, j] == A1A2]
  pheno_A2A2 <- pheno[SNPs[, j] == A2A2]
  
  # Calculate variances for each genotype group
  var_A1A1[j] <- var(pheno_A1A1)
  var_A1A2[j] <- var(pheno_A1A2)
  var_A2A2[j] <- var(pheno_A2A2)
}


# Calculate degrees of freedom for each SNP


df_welch_a = ((var_A1A1/n_A1A1 + var_A2A2/n_A2A2)^2)/((var_A1A1/n_A1A1)^2/(n_A1A1-1) +(var_A2A2/n_A2A2)^2/(n_A2A2-1))
df_welch_d = ((var_A1A2/n_A1A2 + 1/4*(var_A1A1/n_A1A1 + var_A2A2/n_A2A2))^2)/(((var_A1A2/n_A1A2)^2)/(n_A1A2-1)+(var_A1A1/n_A1A1)^2/(4*(n_A1A1-1))+(var_A2A2/n_A2A2)^2/(4*(n_A2A2-1)))

var_a_welch = 0.25*(var_A1A1/n_A1A1 + var_A2A2/n_A2A2)
var_d_welch = var_A1A2/n_A1A2+0.25*(var_A1A1/n_A1A1 + var_A2A2/n_A2A2)

test_a_welch = a_hat/(sqrt(var_a_welch))
test_d_welch = d_hat/(sqrt(var_d_welch))

png(
  filename = file.path(
    out_dir,
    "QQ_plot_test_statistics_a_d.png"
  ),
  width = 1800,
  height = 3200,
  res = 300
)

old_par <- par(
  family = "sans",
  mfrow = c(2, 1),
  mar = c(4.5, 5.5, 3.5, 1),
  oma = c(0, 0, 0, 0)
)

# -------------------------
# Panel A: estimated a
# -------------------------
qqplot(
  qt(
    ppoints(length(test_a_welch)),
    median(df_welch_a)
  ),
  test_a_welch,
  xlab = "Theoretical quantiles",
  ylab = "Sample quantiles",
  pch = 1,
  col = "black",
  lwd = 1.2
)

abline(
  a = 0,
  b = 1,
  col = "red",
  lwd = 1.2
)
mtext(
  "A",
  side = 2,
  line = 4.2,
  at = par("usr")[4],
  las = 1,
  font = 2,
  cex = 1.05
)

# -------------------------
# Panel B: estimated d
# -------------------------
qqplot(
  qt(
    ppoints(length(test_d_welch)),
    median(df_welch_d)
  ),
  test_d_welch,
  xlab = "Theoretical quantiles",
  ylab = "Sample quantiles",
  pch = 1,
  col = "black",
  lwd = 1.2
)

abline(
  a = 0,
  b = 1,
  col = "red",
  lwd = 1.2
)

mtext(
  "B",
  side = 2,
  line = 4.2,
  at = par("usr")[4],
  las = 1,
  font = 2,
  cex = 1.05
)
par(old_par)
dev.off()
