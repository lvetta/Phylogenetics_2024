library(ape)
library(ggtree)
library(ggplot2)
library(plotly)
setwd("C:/Users/Лиза/Desktop/connect")

all_tree <- read.tree("tree_all.newick")
ggtree(all_tree) + 
  geom_text2(aes(label=label), hjust=0, vjust=.15) 
ggsave("all_tree.pdf")
dev.off()

all_tree <- read.tree("tree.newick")

p <- ggtree(all_tree, mrsd = "2013-01-01") + 
  geom_tiplab(size = 2, align = TRUE, linesize = .50) +  
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, vjust = 0.5, size = 3)  


p <- p + theme_tree2() +  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(5, 5, 5, 5, "cm"))  

print(p)


setwd("C:/Users/Лиза/Desktop/connect/Human/Human/")
library(Biostrings)
library(seqinr)
alignment <- read.alignment("Human_all_alg.fasta", format = "fasta")

# Идентификатор для сравнения
reference_id <- "FJ713601.1"  

reference_seq <- NULL

# Поиск референсной последовательности по id
for (i in 1:length(alignment$nam)) {
  if (alignment$nam[i] == reference_id) {
    reference_seq <- alignment$seq[[i]]
    break  
  }
}


# Инициализация вектора для хранения результатов
mutation_counts <- integer(length(alignment$nam) - 1)
names(mutation_counts) <- alignment$nam[alignment$nam != reference_id]  # Имена для остальных последовательностей

# Сравнение с каждой последовательностью
count <- 1
for (i in 1:length(alignment$nam)) {
  if (alignment$nam[i] != reference_id) {
    seq <- alignment$seq[[i]]
    
    
    # Подсчет мутаций
    mutations <- sum(strsplit(reference_seq, NULL)[[1]] != strsplit(seq, NULL)[[1]])
    mutation_counts[count] <- mutations
    count <- count + 1
  }
}

# Вывод результатов
cat("Количество мутаций по сравнению с", reference_id, ":\n")
for (i in seq_along(mutation_counts)) {
  cat(names(mutation_counts)[i], ":", mutation_counts[i], "\n")
}

#определение количества мутаций между неандертальцами, денисевцами и слвременными людьми 
alignment <- read.alignment("old_al.fasta", format = "fasta")

reference_id <- "FJ713601.1"  

reference_seq <- NULL
for (i in 1:length(alignment$nam)) {
  if (alignment$nam[i] == reference_id) {
    reference_seq <- alignment$seq[[i]]
    break  
  }
}
mutation_counts <- integer(length(alignment$nam) - 1)
names(mutation_counts) <- alignment$nam[alignment$nam != reference_id]  # Имена для остальных последовательностей

# Сравнение с каждой последовательностью
count <- 1
for (i in 1:length(alignment$nam)) {
  if (alignment$nam[i] != reference_id) {
    seq <- alignment$seq[[i]]
    
    
    # Подсчет мутаций
    mutations <- sum(strsplit(reference_seq, NULL)[[1]] != strsplit(seq, NULL)[[1]])
    mutation_counts[count] <- mutations
    count <- count + 1
  }
}

# Вывод результатов
cat("Количество мутаций по сравнению с", reference_id, ":\n")
for (i in seq_along(mutation_counts)) {
  cat(names(mutation_counts)[i], ":", mutation_counts[i], "\n")
}
