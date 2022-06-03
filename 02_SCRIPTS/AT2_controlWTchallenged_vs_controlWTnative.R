#####################################################################
#Title: Script Volcano
#Author: Jose A. Ovando
#Date: Marzo/2022

#####################################################################
# Import and select data

library(readr)
library(pillar)
library(dplyr)
library(EnhancedVolcano)

DEG_file <- read_csv("01_DATOS/indmoraV2_CPT1_Cluster1and16and34AT2_controlWTchallenged_vs_controlWTnaive.csv")
colnames(DEG_file)[1] <- "Gene Symbol"
colnames(DEG_file)[3] <- "avg_log2FC"
colnames(DEG_file)[6] <- "adj_P_val"

DEG_file$FC2 <- 2^DEG_file$avg_log2FC
DEG = DEG_file

#####################################################################
#Volcano
#####################################################################

keyvals <- ifelse(DEG_file$avg_log2FC <= -0.2 & DEG_file$adj_P_val <0.05,
                  'blue',  ifelse(DEG_file$avg_log2FC > 0.2 & DEG_file$adj_P_val <0.05,
                                     'red', 'gray'))

#Nombre de valores de cada color
keyvals[is.na(keyvals)] <- "#e3deeb"
names(keyvals)[keyvals == 'red'] <- "Up"
names(keyvals)[keyvals == 'gray'] <- 'Mid'
names(keyvals)[keyvals == 'blue'] <- 'Down'

tiff("volcano_AT2_controlWTchallenged_vs_controlWTnative.tiff", res=300, width = 18, height = 15, units = 'in')
EnhancedVolcano(DEG, x = "avg_log2FC", y = "adj_P_val",
                lab = DEG$`Gene Symbol`,
                xlab = bquote(~avg_log2FC),
                ylab = bquote(~-log[10]~ adj_P_val),
                ylim = c(-0, 150),
                xlim = c(-1.5, 1.5),
                pCutoff = 0.05,
                FCcutoff = 0.2,
                colCustom = keyvals,
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.5,
                axisLabSize = 18,
                pointSize = 4.0,
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabels = c("NS", expression(avg_log2FC), "p-val adj",
                                 expression(adj_P_val ~ and ~ avg_log2FC)),
                legendLabSize = 18,
                title = "",
                titleLabSize = 15,
                subtitleLabSize = 11,
                captionLabSize = 11,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = paste0("N = ", nrow(DEG)),
                subtitle = "",
                labSize = 6,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = T,
                widthConnectors = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#sinnombresdegenes
tiff("volcano_AT2_controlWTchallenged_vs_controlWTnative_sin_nombres.tiff", res=300, width = 18, height = 15, units = 'in')
EnhancedVolcano(DEG, x = "avg_log2FC", y = "adj_P_val",
                lab = DEG$`Gene Symbol`,
                xlab = bquote(~avg_log2FC),
                ylab = bquote(~-log[10]~ adj_P_val),
                ylim = c(-0, 150),
                xlim = c(-1.5, 1.5),
                pCutoff = 0.05,
                FCcutoff = 0.2,
                colCustom = keyvals,
                cutoffLineType = 'solid',
                cutoffLineCol = 'coral4',
                hlineWidth = 0.2,
                cutoffLineWidth = 0.5,
                axisLabSize = 18,
                pointSize = 4.0,
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabels = c("NS", expression(avg_log2FC), "p-val adj",
                                 expression(adj_P_val ~ and ~ avg_log2FC)),
                legendLabSize = 18,
                title = "",
                titleLabSize = 15,
                subtitleLabSize = 11,
                captionLabSize = 11,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                caption = paste0("N = ", nrow(DEG)),
                subtitle = "",
                labSize = 0,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = T,
                widthConnectors = 0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



