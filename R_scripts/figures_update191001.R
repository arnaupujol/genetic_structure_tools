#'#################################################################################
#'#################################################################################
#' Figures of MOI update 191001
#'#################################################################################
#'#################################################################################

## Reads trimming
library(ggplot2)
d <- rbind(AmaMhLDF, CspMhLDF) %>% 
  filter(Dataset == "Set1 Run 2") 

png("QC_seq.png", width = 500, height = 300)
ggplot(d, aes(x = X.BaseNum, y = MisCum*100)) +
   geom_line() +
   scale_x_continuous(name = "Read position", breaks = seq(0, 300, 50)) +
   scale_y_continuous(name = "Cumulate mismatch percentage (%)", breaks = c(0, 0.5, 1, 2, 3, 4, 5, 10)) +
   facet_grid(Marker ~ Type) +
  geom_vline(data = filter(d, Type == "Forward"), aes(xintercept = 275 + 22), colour = "blue", linetype = "dotted") +
  geom_vline(data = filter(d, Type == "Reverse"), aes(xintercept = 253 + 22), colour = "darkblue", linetype = "dotted") +
  theme_bw()
dev.off()

