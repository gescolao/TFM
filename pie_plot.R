#create a datafranme
data <- data.frame(group = c("upregulated", "downregulated"),
                   value = c(279, 79))

#plot the pie plot with the dataframe
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  theme_void() +
  annotate("text", y = 150, x= 1, label= "279 DEGs", 
           color = "white", size = 15) + 
  annotate("text", y = 310, x= 1, label= "79 DEGs", 
           color = "white", size = 15) +
  theme(legend.key.size = unit(2, units = "cm")) + 
  theme(legend.text=element_text(size=15)) + 
  theme(legend.title = element_text(size=0, hjust=0.5))+
  scale_fill_manual(values = c("#430c55","#234F1E"))

                    