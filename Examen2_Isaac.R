#Examen 2: Isaac Flores

#Grafica de volcano
install.packages("packman")
library(pacman)
p_load("readr",
       "ggplot2",
       "matrixTests", 
       "dplyr","ggrepel")

datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/examen2")
head(datos)


controles<- datos %>% 
  filter(Condicion=="Control")
head(controles)

promedio_controles <- controles %>% 
  summarise(Mean_C1=mean(Cx1),
            Mean_C2=mean(Cx2),
            Mean_C3=mean(Cx3),
            Mean_T1=mean(T1),
            Mean_T2=mean(T2),
            Mean_T3=mean(T3) )%>%
  mutate(Gen="Promedio_controles") %>%
  select(7,1,2,3,4,5,6) %>% 
  mutate(Gen="promedio_controles") %>% 
  select(7,1,2,3,4,5,6)
promedio_controles


genes <- datos %>% 
  filter(Condicion == "Target") %>% 
  select(-2)
head(genes)


DCT <- genes %>% 
  mutate(DCT_C1=2^-(Cx1-promedio_controles$Mean_C1),
         DCT_C2=2^-(Cx2-promedio_controles$Mean_C2),
         DCT_C3=2^-(Cx3-promedio_controles$Mean_C3),
         DCT_T1=2^-(T1-promedio_controles$Mean_T1),
         DCT_T2=2^-(T2-promedio_controles$Mean_T2),
         DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)

DCT         

promedio_genes<-DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)

promedio_genes




tvalue_gen <- row_t_welch(promedio_genes[,c("DCT_C1",
                                            "DCT_C2",
                                            "DCT_C3")],
                          promedio_genes[,c("DCT_T1",
                                            "DCT_T2",
                                            "DCT_T3")])
View(tvalue_gen)

FCyPV<- promedio_genes %>% 
  select(1,8,9) %>% 
  mutate(p_value=tvalue_gen$pvalue,
         Fold_change=Mean_DCT_Tx/Mean_DCT_Cx) %>% 
  select(1,4,5)

FCyPV

Logs <- FCyPV %>% 
  mutate(LPV= -log10(p_value),
         LFC= log2 (Fold_change)) %>% 
  select(1,4,5)
Logs


vulcano <- ggplot(Logs, mapping =aes(x=LFC,
                                     y=LPV))+
  geom_point(size=2,
             color= "gray")+
  theme_classic()+
  labs(title="Analisis comparativo de miRNA",
       caption="Isaac Flores",
       x="Log2(2-DDCT",
       y="Log10 (valor de p)")
vulcano  

limite_p <- 0.05

limite_FC <-1.5

down_regulated <- Logs %>% 
  filter(LFC< -log2(limite_FC),
         LPV > -log10(limite_p))
down_regulated

up_regulated <- Logs %>% 
  filter(LFC > log2(limite_FC),
         LPV > -log10(limite_p))

up_regulated

top_down_regulated <- down_regulated %>% 
  arrange(desc(LPV)) %>% 
  head(5)
top_down_regulated

top_up_regulated <- up_regulated %>% 
  arrange(desc(LPV)) %>% 
  head(5)
top_up_regulated





vulcano2 <-vulcano+
  geom_hline(yintercept = -log10(limite_p),linetype="dashed")+
  geom_vline(xintercept = c(-log2(limite_FC), log2(limite_FC)),linetype="dashed")

vulcano2

#derecha sobrepresion

vulcano3 <- vulcano2+
  geom_point(data= up_regulated,
             x = up_regulated$LFC,
             y = up_regulated$LPV,
             alpha= 1,
             size= 3,
             color="#E64835B2")+
  geom_point(data= down_regulated,
             x = down_regulated$LFC,
             y = down_regulated$LPV,
             alpha= 1,
             size= 3,
             color="#135477B2")

vulcano3

vulcano4 <- vulcano3+geom_label_repel (data = top_up_regulated,
                                       mapping = aes(x = top_up_regulated$LFC,
                                                     y = top_up_regulated$LPV),
                                       label = top_up_regulated$miRNA,
                                       max.overlaps=50)+
  geom_label_repel(data = top_down_regulated,
                   mapping = aes(x = top_down_regulated$LFC,
                                 y = top_down_regulated$LPV),
                   label = top_down_regulated$miRNA,
                   max.overlaps=50)

vulcano4

ggsave("Examen2Volcano.jpeg",
       plot = vulcano4,
       height = 5,width = 6,
       dpi = 300)

