# Load the necessary package
library(VennDiagram)

# Read the data from files
data1 <- read.table("DEGs94.txt")
data2 <- read.table("DEGs02.txt")
data3 <- read.table("DEGs23.txt")
data4 <- read.table("DEGs49.txt")
data5 <- read.table("DEGs71.txt")
data6 <- read.table("DEGs12.txt")


degs.res_UP1=data1[data1$logFC >= 1   & data1$adj.P.Val <= 0.05,]  
degs.res_DOWN1=data1[data1$logFC  <= -1   & data1$adj.P.Val <= 0.05,] 

degs.res_UP2=data2[data2$logFC >= 1   & data2$adj.P.Val <= 0.05,]  
degs.res_DOWN2=data2[data2$logFC  <= -1   & data2$adj.P.Val <= 0.05,] 

degs.res_UP3=data3[data3$logFC >= 1   & data3$adj.P.Val <= 0.05,]  
degs.res_DOWN3=data3[data3$logFC  <= -1   & data3$adj.P.Val <= 0.05,] 

degs.res_UP4=data4[data4$logFC >= 1   & data4$adj.P.Val <= 0.05,]  
degs.res_DOWN4=data4[data4$logFC  <= -1   & data4$adj.P.Val <= 0.05,] 

degs.res_UP5=data5[data5$logFC >= 1   & data5$adj.P.Val <= 0.05,]  
degs.res_DOWN5=data5[data5$logFC  <= -1   & data5$adj.P.Val <= 0.05,] 

degs.res_UP6=data6[data6$logFC >= 1   & data6$adj.P.Val <= 0.05,]  
degs.res_DOWN6=data6[data6$logFC  <= -1   & data6$adj.P.Val <= 0.05,] 


# Extract the DEG names
up_names   <- c(rownames(degs.res_UP1),rownames(degs.res_UP2),rownames(degs.res_UP3)
                ,rownames(degs.res_UP4),rownames(degs.res_UP5),rownames(degs.res_UP6) )
down_names <- c(rownames(degs.res_DOWN1),rownames(degs.res_DOWN2),rownames(degs.res_DOWN3)
                ,rownames(degs.res_DOWN4),rownames(degs.res_DOWN5),rownames(degs.res_DOWN6) )

# Find the counts of each
up_counts <- table(up_names)
down_counts <- table(down_names)


# Find the names that appear in at least two files
up_overlap <- up_counts[up_counts >= 2]
down_overlap <- down_counts[down_counts >= 2]

write.table(up_overlap, "Overlapped_up_DEGs.txt")
write.table(down_overlap, "Overlapped_down_DEGs.txt")



#[2] Visualization
# [A] ----------- Volcano (UP/Down)-----------  

All <- dplyr::bind_rows(data1, data2, data3, data4, data5, data6)

res2 <- All %>%
  mutate(gene_type = case_when(logFC >= 1  & adj.P.Val <= 0.05 ~ "up",
                               logFC <= -1  & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
cols <- c("up" = "#FF0000", "down" = "#00FF00", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

res2 %>%
  ggplot(aes(x = (logFC),
             y = -log10(adj.P.Val),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = (0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10)) 

#[2] Visualization




