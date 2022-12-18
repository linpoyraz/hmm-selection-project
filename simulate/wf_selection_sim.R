# simulates haploid wright fisher population
# sample 100 haploid individuals every 10 generations 
# usage: wf_selection_sim.R <selection_file> <Ne> <initial frequency> 
library(data.table)

set.seed(1)

args = commandArgs(trailingOnly = T)
Ne = as.numeric(args[2])
init_freq = as.numeric(args[3])

# read selection file: first column will be time (generation), and second selection coef
selection_coef = read.table(args[1], header = T)
names(selection_coef) = c("t", "s")

# initialize results data frame. initialize true AF data frame 
result = data.frame("t" = NaN, "n" = NaN, "a" = NaN, "init_freq" = init_freq, "Ne" = Ne, "s" = NaN, "freq" = NaN)
true = data.frame("t" = NaN, "freq" = NaN, "init_freq" = init_freq, "Ne" = Ne, "s" = NaN)
# begin WF simulation: the number of generations simulated will be equal to the number
# specified in the selection coefficient file 
freq = init_freq
for (i in 1:nrow(selection_coef)){
  s = selection_coef$s[i]
  formula = ((1 + s) * freq) / (1 + (s * freq))
  freq = sum(rbinom(n = Ne, size = 1, prob = formula)) / Ne
  gen_true = data.frame("t" = i, "freq" = freq, "init_freq" = init_freq, "Ne" = Ne, "s" = s)
  true = rbind(true, gen_true)
  sample = sum(rbinom(n = 100, size = 1, prob = freq))
  gen_result = data.frame("t" = i, "n" = ifelse(i %% 10 == 0 | i == 1, 100, 0), "a" = ifelse(i %% 10 == 0 | i == 1, sample, 0), "init_freq" = init_freq, "Ne" = Ne, "s" = s, "freq" = freq)
  result = rbind(result, gen_result)
}

# remove first fake row 
result = result[2:nrow(result),]
true = true[2:nrow(true),]
# write file to simulated_files directory 
f = paste("./simulated_files/Ne_", Ne, "_init_freq_", init_freq, "_", args[1], sep = "")
fwrite(result, f, quote = F, sep = "\t", row.names = F, col.names = T)

f = paste("./simulated_files/true_Ne_", Ne, "_init_freq_", init_freq, "_", args[1], sep = "")
fwrite(true, f, quote = F, sep = "\t", row.names = F, col.names = T)





