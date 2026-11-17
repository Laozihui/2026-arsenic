library(dplyr)
library(readxl)
library(stringr)

data <- read.csv("14n.csv", header=TRUE,check.names = F)
library(mediation)
library(parallel);#加载并行计算包

cl <- makeCluster(10);# 初始化cpu集群

clusterEvalQ(cl,library(mediation));#添加并行计算中用到的包

clusterExport(cl,'data');

do.data= function(.x){
  lapply(colnames(data)|>stringr::str_extract("^XQ_.*")|>na.omit(), function(.y){
    lapply(colnames(data)|>stringr::str_extract("^g_.*")|>na.omit(), function(.z){
      
      data[,c(.z,.x,.y)]->datass
      
      colnames(datass)=c("X","Y","M")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      set.seed(1)
      summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      colnames(datass)=c("X","M","Y")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      summary2=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      data.frame(group=paste0(.z,"==>",.y,"==>",.x),
                 Pval_mediate=summary$d.avg.p,
                 Pval_direct=summary$z.avg.p,
                 Pval_mediate_inverse=summary2$d.avg.p,
                 Pval_direct_inverse=summary2$z.avg.p
      )
      #rownames(k)=paste0(.z,"==>",.y,"==>",.x,"===", rownames(k))
    })
  })
}
h=colnames(data)%>%str_extract("^ZZ_.*")%>%na.omit
get1=parLapply(cl,h , do.data)
zfs_get=lapply(get1,\(x) lapply(x,\(.x) {.x%>%do.call(rbind,.) })%>%{do.call(rbind,.) })%>%do.call(rbind,.)
zfs_get$Qval_mediate=p.adjust(zfs_get$Pval_mediate,method = "BH")
zfs_get$Qval_mediate_inverse=p.adjust(zfs_get$Pval_mediate_inverse,method = "BH")
zfs_get=zfs_get[which(zfs_get$Qval_mediate<0.05&zfs_get$Pval_mediate_inverse>0.05),]

zfs_get$input= zfs_get$group%>%str_extract("^.*?(?===>)")
zfs_get$line= zfs_get$group%>%str_extract("(?<=>).*?(?===>)")
zfs_get$output= zfs_get$group%>%str_extract("(?<=>)ZZ_.*")
net=zfs_get
net$fre=1
write.csv(net, "net14.csv", row.names = TRUE)
library(ggplot2)
library(ggalluvial)
##install.packages("ggalluvial")
library(gridExtra)

p<-ggplot(net,aes(axis1 = net$input, axis2 = net$line, axis3 = net$output,y= net$fre))+
  scale_x_discrete(limits = c("Microbiome", "XQ", "zz")) +
  geom_alluvium(aes(fill = net$line),cex=0.1)+
  geom_stratum(width=0.4) +
  
  
  theme_minimal()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  
  theme(legend.position="none",axis.text = element_text(size = 10),
        axis.title = element_text(size = 6),panel.grid=element_blank())
  #scale_fill_manual(values = c("#6DAFD1", "#B72230", "#FFC857", "#FF8C42", "#76B041", "#996633", "#845EC2"))#+
p

net<-read.csv("net14.csv",header = T,check.names = F)


#data1<-t(data[,-1])
#write.csv(data1, "new_data14.csv", row.names = TRUE)
#r<-read.csv("new_data14.csv",header = F,check.names = F)

net_plot <- function(net) {
  
  nodes <- as.character(
    unlist(
      net[c("input", "line", "output")],
      use.names = FALSE
    )
  )
  
  df <- data.frame(
    x = c(0, 0.5, 1),
    y = c(0, 1, 0),
    group = nodes
  )
  
  pairs <- combn(nodes, 2)
  
  f <- apply(
    pairs,
    2,
    function(x) {
      if (all(x %in% colnames(data))) {
        k <- cor.test(
          data[[x[1]]],
          data[[x[2]]]
        )
        
        data.frame(
          come = x[1],
          to = x[2],
          cor = unname(k$estimate),
          p = k$p.value
        )
      } else {
        NULL
      }
    }
  ) %>%
    do.call(rbind, .)
  
  f1 <- f[
    f$come %in% df$group[1:2] &
      f$to %in% df$group[1:2],
  ]
  
  f2 <- f[
    f$come %in% df$group[2:3] &
      f$to %in% df$group[2:3],
  ]
  
  f3 <- f[
    f$come %in% df$group[c(1, 3)] &
      f$to %in% df$group[c(1, 3)],
  ]
  
  if (any(c(f1$p, f2$p, f3$p) < 0.05)) {
    
    line <- tibble(
      x = c(0, 0.5, 1, 0.0001),
      y = c(0, 1, 0, 0.00001)
    )
    
    ggplot() +
      geom_path(
        data = line,
        aes(x = x, y = y),
        linewidth = 10,
        colour = "grey"
      ) +
      ggplot2::annotate(
        x = 0.1,
        y = 0.5,
        angle = 65,
        label = paste0(
          "cor=", round(f1$cor, 4),
          " , pval=", round(f1$p, 4)
        ),
        geom = "text"
      ) +
      ggplot2::annotate(
        x = 0.1 + 0.5 + 0.2,
        y = 0.6,
        angle = 295,
        label = paste0(
          "cor=", round(f2$cor, 4),
          " , pval=", round(f2$p, 4)
        ),
        geom = "text"
      ) +
      ggplot2::annotate(
        x = 0.5,
        y = -0.2,
        angle = 0,
        label = paste0(
          "cor=", round(f3$cor, 4),
          " , pval=", round(f3$p, 4)
        ),
        geom = "text"
      ) +
      geom_point(
        data = df,
        aes(x = x, y = y),
        size = 10
      ) +
      coord_fixed() +
      ggrepel::geom_text_repel(
        data = df[1, ],
        aes(x = x, y = y - 0.1, label = group)
      ) +
      ggrepel::geom_text_repel(
        data = df[2, ],
        aes(x = x, y = y + 0.1, label = group)
      ) +
      ggrepel::geom_text_repel(
        data = df[3, ],
        aes(x = x, y = y - 0.1, label = group)
      ) +
      theme_void()
  } else {
    NULL
  }
}

plot_data=1:nrow(net)%>%lapply(\(.x){
  net_plot(net[.x,])}
)
dir.create("plots14")

# 遍历每个绘图并保存
for (i in seq_along(plot_data)) {
  ggsave(paste0("plots14/plot_", i, ".pdf"), plot = plot_data[[i]], width = 8, height = 6, units = "in")
}
