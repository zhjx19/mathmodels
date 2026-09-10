devtools::load_all("E:/Github/mathmodels")
data3 = data.frame(industry=c('A','A','A','A','B','B','B','B'),
  area=c('East','East','West','West','East','East','West','West'),
  province=c('Shanghai','Beijing','Sichuan','Yunnan','Shanghai','Beijing','Sichuan','Yunnan'),
  avg_wage=c(500,400,80,60,300,250,50,40), emp_num=c(100,80,90,70,120,100,80,60))
d = data3; names(d) = c('g1','g2','prov','y','pop')
Ttot = theil(d$y, d$pop)
g1agg = d |> dplyr::summarise(y=sum(y*pop)/sum(pop), pop=sum(pop), .by=g1) |>
  dplyr::mutate(s=y*pop/sum(y*pop), p=pop/sum(pop))
Tb_g1 = theil(g1agg$y, g1agg$pop)
res = NULL
for(i in seq_len(nrow(g1agg))){
  g = g1agg$g1[i]; di = d[d$g1==g,]
  a2 = di |> dplyr::summarise(y=sum(y*pop)/sum(pop), pop=sum(pop), .by=g2) |>
    dplyr::mutate(s=y*pop/sum(y*pop), p=pop/sum(pop))
  Tb_g2 = sum(a2$s*log(a2$s/a2$p))
  Tw_g2 = 0
  for(j in seq_len(nrow(a2))){
    gg = a2$g2[j]; dc = di[di$g2==gg,]
    sc = sum(dc$y*dc$pop)/sum(di$y*di$pop)
    Tw_g2 = Tw_g2 + sc*theil(dc$y, dc$pop)
  }
  res = rbind(res, data.frame(g=g, s=g1agg$s[i], Tb_g2=Tb_g2, Tw_g2=Tw_g2, T_g1=Tb_g2+Tw_g2))
}
cat("Ttot direct         =", format(Ttot, digits=8), "\n")
cat("Tb_g1               =", format(Tb_g1, digits=8), "\n")
cat("sum s_g1 * T_g1     =", format(sum(res$s*res$T_g1), digits=8), "\n")
cat("reconstructed total =", format(Tb_g1+sum(res$s*res$T_g1), digits=8), "\n")
cat("package TT          =", format(theil_g2_cross(data3,'industry','area','avg_wage','emp_num')$theil['T'], digits=8), "\n")
print(res)
