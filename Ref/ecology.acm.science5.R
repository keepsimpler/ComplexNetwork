###############################################################################
####  f + g are fixed
###############################################################################
# Explore the equations
sigma = 0.05  # sd of half normal distribution
EX = sigma * sqrt(2/pi)  # mean of half normal distribution E(|X|)
EX2 = sigma^2  # E(|X|^2)
E <- expression(C * (f * Pm - g * (1 - Pa - Pm)) * EX)  # mean of off-diagonal elements
# second moment of off-diagonal elements
E2UP <- expression(C * (f^2 * Pm + h^2 * Pa + g^2 * (1 - Pa - Pm)) * EX2)
# mean of the product of interaction pairs
E2DOWN <- expression(C * (f^2 * Pm - h^2 * Pa + g^2 * (1 - Pa - Pm)) * EX^2)
V <- parse(text = paste(E2UP, '- (', E, ')^2'))  # E2UP - E^2
VC <- parse(text = paste(E2DOWN, '- (', E, ')^2'))  # E2DOWN - E^2
RHO <- parse(text = paste('(', VC, ') / (', V, ')'))  # VC / V
# dot eigenvalue: (S - 1) * E - s
Rs <- parse(text = paste('(S - 1) *', E, ' - s'))
# ellipse eigenvalue: - E - s + sqrt(S * V) * (1 + rho)
Re <- parse(text = paste('-', E, '- s', '+ sqrt(S * (', V, ') ) * (1 + (', RHO, ') )'))
Re2 <- parse(text = paste('-', E, '- s', '+ sqrt(S * (', V, ') ) + sqrt(S/(', V, ') ) * (', VC, ')'))  # - E - s + sqrt(S * V) + (sqrt(S/V) * VC)
Re1deriv <- D(Re, 'Pm')  # first derivative of Re
Re1deriv2 <- D(Re2, 'Pm')  # first derivative of Re
Re.sub.Rs <- parse(text = paste(Re, '- (', Rs, ')'))
TS <- 1

# define the environment used by eval function
f_env <- function(Pa, Pm, Pc, C, alpha) {
  #amplification = TS / (Pm^(1+alpha) + Pc^(1+alpha) + Pa^(1+alpha))
  #f = amplification * Pm^alpha
  #g = amplification * Pc^alpha
  #h = amplification * Pa^alpha
  if (round(Pm, 5) == 0) f = 0
  else f = Pm^alpha
  if (round(Pc, 5) == 0) g = 0
  else g = Pc^alpha
  if (round(Pa, 5) == 0) h = 0
  else h = Pa^alpha
  amplification = TS / (f*Pm + g*Pc + h*Pa)
  f = amplification * f
  g = amplification * g
  h = amplification * h
  env <- list(
    Pa = Pa,
    C = C,
    f = f,
    g = g,
    h = h,
    EX = EX,
    EX2 = EX2,
    Pm = Pm,
    Pc = Pc
  )
  return(env)
}

# coefficients
S = 100  # species number
s = 1  # self-regulation
C = seq(from = 0.1, to = 0.9, by = 0.4)  # connectance
Pm = seq(from = 0.0, to = 1., by = 0.01)  # antagonism proportion
Pc = seq(from = 0.0, to = 1., by = 0.01)  #  proportion
Pa = seq(from = 0.0, to = 1., by = 0.01)  #  proportion
coeffs <- expand.grid(C = C, Pm = Pm, Pc = Pc, Pa = Pa)
coeffs <- subset(coeffs, round(Pm + Pc + Pa, 5) == 1)
alpha = data.frame(alpha = seq(from = -1., to = 0., by = 0.1)) # c(-1, -0.5, 0, 0.5, 1)
coeffs = merge(coeffs, alpha)
coeffs$id = 1:nrow(coeffs)


##### Simulate to get the measurements [Re, Rs, R, D(Re)/D(Pm), Re-Rs]
##### according to (Pa,Pc,Pm) and (C,f)
library(plyr)
ret.h0.v5 <- ddply(coeffs, .variables = .(id), function(coeff) {
  with(coeff, {
    print(id)
    if (round(Pm, 5) == 0) f = 0
    else f = Pm^alpha
    if (round(Pc, 5) == 0) g = 0
    else g = Pc^alpha
    if (round(Pa, 5) == 0) h = 0
    else h = Pa^alpha
    amplification = TS / (f*Pm + g*Pc + h*Pa)
    f = amplification * f
    g = amplification * g
    h = amplification * h
    env <- f_env(Pa, Pm, Pc, C, alpha)
    #mean <- eval(E, envir = env)
    #mean2up <- eval(E2UP, envir = env)
    #mean2down <- eval(E2DOWN, envir = env)
    #v <- eval(V, envir = env)
    #vc <- eval(VC, envir = env)
    #rho <- eval(RHO, envir = env)
    re <- eval(Re, envir = env)
    rs <- eval(Rs, envir = env)
    #re.sub.rs <- eval(Re.sub.Rs, envir = env)
    lev <- pmax(re, rs)
    #re1deriv <- eval(Re1deriv, envir = env)
    #re1deriv2 <- eval(Re1deriv2, envir = env)
    #cbind(Pa, C, alpha, f, g, h, Pm, Pc, mean, mean2up, mean2down, v, vc, rho, re, rs, lev, re.sub.rs, re1deriv, re1deriv2)
    cbind(Pa, Pm, Pc, C, alpha, f, g, h, amplification, re, rs, lev)
  })
})
all(round(ret.h0.v5$Pm * ret.h0.v5$f + ret.h0.v5$Pc * ret.h0.v5$g + ret.h0.v5$Pa * ret.h0.v5$h,5) == TS)

# Fig1: mutualistic-competitive
tmp <- subset(ret.h0.v5, C == 0.5 & round(Pa,5) == 0.0 & round(Pm,5) > 0 &  round(Pc,5) > 0 & round(alpha, 5) %in% c(-1, -0.8, -0.6, -0.4, -0.2, 0))
tmp <- tmp[, c('Pa', 'Pm', 'Pc', 'alpha', 'f', 'g', 'h', 'lev')]
library(data.table)
tmp2 <- data.table(tmp)
tmp2 <- tmp2[, .SD[which.min(lev)], by = alpha]
p.lev <- ggplot(data = tmp, mapping = aes(x = Pm, y = lev, colour = factor(alpha))) +  # , group = factor(f)
  geom_line() +
  geom_point(data = tmp2, aes(x = Pm, y = lev, group = factor(alpha)), size = 2, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  guides(colour = guide_legend(title = expression(alpha), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
#  theme(legend.position = c(.5, .5)) +
  theme_bw()
p.lev

# Fig1: antagonistic-mutualistic
tmp <- subset(ret.h0.v5, C == 0.5 & round(Pc,5) == 0.0 & round(Pm,5) > 0 &  round(Pa,5) > 0 & round(alpha, 5) %in% c(-1, -0.8, -0.6, -0.4, -0.2, 0)) #
tmp <- tmp[, c('Pa', 'Pm', 'Pc', 'alpha', 'f', 'g', 'h', 'lev')]
library(data.table)
tmp2 <- data.table(tmp)
tmp2 <- tmp2[, .SD[which.min(lev)], by = alpha]
p.lev <- ggplot(data = tmp, mapping = aes(x = Pa, y = lev, colour = factor(alpha))) +  # , group = factor(f)
  geom_line() +
  geom_point(data = tmp2, aes(x = Pa, y = lev, group = factor(alpha)), size = 2, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[a]), y = expression(R)) +
  guides(colour = guide_legend(title = expression(alpha), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev

# Fig1: competitive-antagonistic
tmp <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) == 0.0 & round(Pc,5) > 0 &  round(Pa,5) > 0 & round(alpha, 5) %in% c(-1, -0.8, -0.6, -0.4, -0.2, 0)) #
tmp <- tmp[, c('Pa', 'Pm', 'Pc', 'alpha', 'f', 'g', 'h', 'lev')]
library(data.table)
tmp2 <- data.table(tmp)
tmp2 <- tmp2[, .SD[which.min(lev)], by = alpha]
p.lev <- ggplot(data = tmp, mapping = aes(x = Pa, y = lev, colour = factor(alpha))) +  # , group = factor(f)
  geom_line() +
  geom_point(data = tmp2, aes(x = Pa, y = lev, group = factor(alpha)), size = 2, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[a]), y = expression(R)) +
  guides(colour = guide_legend(title = expression(alpha), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev

#Fig.2
# Fig. 3
c05alpha10 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -1 ) #
write.table(c05alpha10[,c('Pm','Pa','Pc','lev')], file = 'c05alpha10.csv', sep = ',',row.names = F, col.names = F)
c05alpha08 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -0.8)
write.table(c05alpha08[,c('Pm','Pa','Pc','lev')], file = 'c05alpha08.csv', sep = ',',row.names = F, col.names = F)
c05alpha06 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -0.6)
write.table(c05alpha06[,c('Pm','Pa','Pc','lev')], file = 'c05alpha06.csv', sep = ',',row.names = F, col.names = F)
c05alpha04 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -0.4)
write.table(c05alpha04[,c('Pm','Pa','Pc','lev')], file = 'c05alpha04.csv', sep = ',',row.names = F, col.names = F)
c05alpha02 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -0.2)
write.table(c05alpha02[,c('Pm','Pa','Pc','lev')], file = 'c05alpha02.csv', sep = ',',row.names = F, col.names = F)
c05alpha00 <- subset(ret.h0.v5, C == 0.5 & round(Pm,5) > 0 & round(Pa,5) >= 0 & round(Pc,5) > 0 & round(alpha, 5) == -0.0)
write.table(c05alpha00[,c('Pm','Pa','Pc','lev')], file = 'c05alpha00.csv', sep = ',',row.names = F, col.names = F)

c05alpha10[c05alpha10$lev == min(c05alpha10$lev), ]
c05alpha08[c05alpha08$lev == min(c05alpha08$lev), ]
c05alpha06[c05alpha06$lev == min(c05alpha06$lev), ]

c01alpha10 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -1)
write.table(c01alpha10[,c('Pm','Pa','Pc','lev')], file = 'c01alpha10.csv', sep = ',',row.names = F, col.names = F)
c01alpha08 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -0.8)
write.table(c01alpha08[,c('Pm','Pa','Pc','lev')], file = 'c01alpha08.csv', sep = ',',row.names = F, col.names = F)
c01alpha06 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -0.6)
write.table(c01alpha06[,c('Pm','Pa','Pc','lev')], file = 'c01alpha06.csv', sep = ',',row.names = F, col.names = F)
c01alpha04 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -0.4)
write.table(c01alpha04[,c('Pm','Pa','Pc','lev')], file = 'c01alpha04.csv', sep = ',',row.names = F, col.names = F)
c01alpha02 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -0.2)
write.table(c01alpha02[,c('Pm','Pa','Pc','lev')], file = 'c01alpha02.csv', sep = ',',row.names = F, col.names = F)
c01alpha00 <- subset(ret.h0.v5, C == 0.1 & round(alpha, 5) == -0.0)
write.table(c01alpha00[,c('Pm','Pa','Pc','lev')], file = 'c01alpha00.csv', sep = ',',row.names = F, col.names = F)

c09alpha10 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -1)
write.table(c09alpha10[,c('Pm','Pa','Pc','lev')], file = 'c09alpha10.csv', sep = ',',row.names = F, col.names = F)
c09alpha08 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -0.8)
write.table(c09alpha08[,c('Pm','Pa','Pc','lev')], file = 'c09alpha08.csv', sep = ',',row.names = F, col.names = F)
c09alpha06 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -0.6)
write.table(c09alpha06[,c('Pm','Pa','Pc','lev')], file = 'c09alpha06.csv', sep = ',',row.names = F, col.names = F)
c09alpha04 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -0.4)
write.table(c09alpha04[,c('Pm','Pa','Pc','lev')], file = 'c09alpha04.csv', sep = ',',row.names = F, col.names = F)
c09alpha02 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -0.2)
write.table(c09alpha02[,c('Pm','Pa','Pc','lev')], file = 'c09alpha02.csv', sep = ',',row.names = F, col.names = F)
c09alpha00 <- subset(ret.h0.v5, C == 0.9 & round(alpha, 5) == -0.0)
write.table(c09alpha00[,c('Pm','Pa','Pc','lev')], file = 'c09alpha00.csv', sep = ',',row.names = F, col.names = F)


tmp <- subset(ret.h0.v5, C == 0.5  & alpha == -0.5 & round(Pm, 5) > 0 & round(Pc, 5) > 0 & (round(Pa, 5) == 0 | round(Pa, 5) == 0.01))
ggplot(tmp, aes(x = Pm, y = lev, group = factor(Pa), color = factor(Pa))) + geom_line()

tmp <- subset(ret.h0.v5, C == 0.5  & alpha == -1 & round(Pa, 5) > 0 & round(Pm, 5) > 0 & (round(Pc, 5) == 0 | round(Pc, 5) == 0.01))
ggplot(tmp, aes(x = Pa, y = lev, group = factor(Pc), color = factor(Pc))) + geom_line()

tmp <- subset(ret.h0.v5, C == 0.5  & alpha == -0.5 & round(Pa, 5) > 0 & round(Pc, 5) > 0 & (round(Pm, 5) == 0 | round(Pm, 5) == 0.01))
ggplot(tmp, aes(x = Pa, y = lev, group = factor(Pm), color = factor(Pm))) + geom_line()


tmp <- subset(ret.h0.v5, C == 0.5  & alpha == -0.5)
tmp <- tmp[, c('Pa', 'Pm', 'Pc', 'f', 'g', 'h', 'lev')]
plot.ggtern.2(tmp,'lev')

plot.ggtern.2 <- function(df, colname) {
  require(ggtern)
  # sequence: Pa --> Pc --> Pm
  ggtern(df,aes(Pa, Pc, Pm)) +
    geom_point(aes_string(fill=colname),color="black",shape=23,size=3) +
    scale_fill_gradient(low="blue",high="red") +
    theme(legend.position=c(0,1),legend.justification=c(0,1)) +
    labs(fill=colname, x = expression(P[a]), y = expression(P[c]), z = expression(P[m])) +
    theme_tern_bw()
}

###############################################################################

f_Re1deriv <- function(Pm, Pa, C, f, g, h) {
  env <- list(
    Pa = Pa,
    C = C,
    f = f,
    g = g,
    h = h,
    EX = EX,
    EX2 = EX2,
    Pm = Pm
  )
  eval(Re1deriv, envir = env)
}
f_Re.sub.Rs <- function(Pm, Pa, C, f, g, h) {
  env <- list(
    Pa = Pa,
    C = C,
    f = f,
    g = g,
    h = h,
    EX = EX,
    EX2 = EX2,
    Pm = Pm
  )
  eval(Re.sub.Rs, envir = env)
}

S = 100  # species number
s = 1  # self-regulation
Pa = seq(from = 0, to = 0.9, by = 0.1)  # antagonism proportion
C = seq(from = 0.1, to = 1, by = 0.01)  # connectance
epsilon = seq(from = 0.1, to = 1.5, by = 0.01)  # e = f / g
h = 1 # seq(from = 0., to = 1., by = 0.1)  # antagonism interaction strength
coeffs <- expand.grid(Pa = Pa, C = C, epsilon = epsilon, h = h)
coeffs$f <- FG * coeffs$epsilon / (coeffs$epsilon + 1)  # f = FG * e / (e + 1)
coeffs$g <- FG / (coeffs$epsilon + 1)  # g = FG / (e + 1)
#f = seq(from = 0., to = 1., by = 0.1)  # cooperation interaction strength
#g = seq(from = 0., to = 1., by = 0.1)  # competition interaction strength
#h = seq(from = 0., to = 1., by = 0.1)  # antagonism interaction strength
#coeffs <- expand.grid(Pa = Pa, C = C, f = f, g = g, h = h)
#coeffs <- subset(coeffs, round(f + g + h, 5) == 1)
#coeffs <- subset(coeffs, round(f, 5) != 1 & round(g, 5) != 1 & round(h, 5) != 1)
coeffs$id = 1:nrow(coeffs)

ret.h0.diff <- ddply(coeffs, .variables = .(id), function(coeff) {
  with(coeff, {
    print(id)
    if (f_Re1deriv(0, Pa, C, f, g, h) > 0 &
        f_Re1deriv(1 - Pa, Pa, C, f, g, h) < 0) {  # when Pm=0, first derivative must be positive
      re.max <- uniroot(f_Re1deriv, Pa, C, f, g, h, interval=c(0, 1 - Pa), tol = 1e-5)$root
    } else if (f_Re1deriv(0, Pa, C, f, g, h) < 0) {
      re.max <- 0
    } else {
      re.max <- 1 - Pa
    }
    if (f_Re.sub.Rs(0, Pa, C, f, g, h) > 0 &
        f_Re.sub.Rs(1 - Pa, Pa, C, f, g, h) < 0) {  # when Pm=0, first derivative must be positive
      re.equal.rs <- uniroot(f_Re.sub.Rs, Pa, C, f, g, h, interval=c(0, 1 - Pa), tol = 1e-5)$root
    } else if (f_Re.sub.Rs(0, Pa, C, f, g, h) < 0) {
      re.equal.rs <- 0
    }
    else {
      re.equal.rs <- 1 - Pa
    }
    cbind(Pa, C, epsilon, f, g, h, re.max, re.equal.rs)
  })
})
ret.h0.diff$diff <- (ret.h0.diff$re.equal.rs - ret.h0.diff$re.max)/(1-ret.h0.diff$Pa)

f_Re <- function(Pm, Pa, C, f, g, h) {
  env <- list(
    Pa = Pa,
    C = C,
    f = f,
    g = g,
    h = h,
    EX = EX,
    EX2 = EX2,
    Pm = Pm
  )
  eval(Re, envir = env)
}
f_Rs <- function(Pm, Pa, C, f, g, h) {
  env <- list(
    Pa = Pa,
    C = C,
    f = f,
    g = g,
    h = h,
    EX = EX,
    EX2 = EX2,
    Pm = Pm
  )
  eval(Rs, envir = env)
}

ret.h0.diff$re.max.value <- with(ret.h0.diff, {
  Pm <- re.max
  f_Re(Pm, Pa, C, f, g, h)
})
ret.h0.diff$re.equal.rs.value <- with(ret.h0.diff, {
  Pm <- re.equal.rs
  f_Re(Pm, Pa, C, f, g, h)
})
ret.h0.diff$rs.max.value <- with(ret.h0.diff, {
  Pm <- 1 - Pa
  f_Rs(Pm, Pa, C, f, g, h)
})
ret.h0.diff$re.pm0.value <- with(ret.h0.diff, {
  Pm <- 0
  f_Re(Pm, Pa, C, f, g, h)
})

ret.h0.diff$deep1 <- (ret.h0.diff$re.max.value - ret.h0.diff$re.equal.rs.value)/(ret.h0.diff$re.max.value - ret.h0.diff$re.pm0.value)
ret.h0.diff[ret.h0.diff$diff <= 0, "deep1"] <- 0

ret.h0.diff$deep <- (ret.h0.diff$re.max.value - ret.h0.diff$re.equal.rs.value)/(ret.h0.diff$rs.max.value - ret.h0.diff$re.equal.rs.value)
#subset(ret.h0.diff, rs.max.value < re.equal.rs.value & diff >= 0 & deep >= 0)
# when the interval vanishs, the deep of interval will also vanish
ret.h0.diff[ret.h0.diff$diff <= 0, "deep"] <- 0
# several exceptional cases
# 1, diff == 0, re.max == re.equal.rs == 1 - Pa, deep == 0
#subset(ret.h0.diff, rs.max.value < re.equal.rs.value & deep == 0)
# 2,  rs.max.value < re.equal.rs.value & deep < 0, should be replaced by 1e10
#subset(ret.h0.diff, rs.max.value < re.equal.rs.value & deep < 0 & re.max.value > re.equal.rs.value)
ret.h0.diff[ret.h0.diff$rs.max.value < ret.h0.diff$re.equal.rs.value & ret.h0.diff$deep < 0 & ret.h0.diff$re.max.value > ret.h0.diff$re.equal.rs.value, "deep"] <- 1e10

# average of D1 + D3
#a fake rs.max.value, in order to facilite calculation of D3 when rs.max.value < rs.equal.rs.value
ret.h0.diff$rs.max.value.fake <- ret.h0.diff$rs.max.value
ret.h0.diff[ret.h0.diff$rs.max.value - ret.h0.diff$re.equal.rs.value < 0 , "rs.max.value.fake"] = ret.h0.diff[ret.h0.diff$rs.max.value - ret.h0.diff$re.equal.rs.value < 0 , "re.equal.rs.value"]

ret.h0.diff$deepavg <- (ret.h0.diff$re.max.value - ret.h0.diff$re.equal.rs.value)/((ret.h0.diff$re.max.value - ret.h0.diff$re.pm0.value + ret.h0.diff$rs.max.value.fake - ret.h0.diff$re.equal.rs.value)/2)
ret.h0.diff[ret.h0.diff$diff <= 0, "deepavg"] <- 0


ret.h0.diff.pa0 <- subset(ret.h0.diff, Pa == 0 )

# Figure 1:
# create a figure to show the trend of [re] and [rs] as a function of [Pm]
# such a function are affected by [Pa], [C], [f]
ret.h0.pa08.c08.h10 <- subset(ret.h0, Pa == 0. & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1)
ret.h0.pa08.c08.h10 <- subset(ret.h0.pa08.c08.h10, round(epsilon, digits = 5) %in% round(seq(0.1,1.5,0.1),digits=5))
ret.h0.diff.pa08.c08.h10 <- subset(ret.h0.diff, Pa == 0. & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1)
ret.h0.diff.pa08.c08.h10 <- subset(ret.h0.diff.pa08.c08.h10, round(epsilon, digits = 5) %in% round(seq(0.1,1.5,0.1),digits=5))

# Fig. 1
re.min <- min(ret.h0.pa08.c08.h10$re)
re.max <- max(ret.h0.pa08.c08.h10$re)
p.lev <- ggplot(data = ret.h0.pa08.c08.h10, mapping = aes(x = Pm, y = lev, colour = factor(epsilon))) +  # , group = factor(f)
  geom_line(aes(x = Pm, y = re, colour = factor(epsilon)), linetype = 3) +
  geom_line() +
  geom_point(data = ret.h0.diff.pa08.c08.h10, aes(x = re.equal.rs, y = re.equal.rs.value, group = factor(epsilon)), size = 2, shape = 23, colour = 'black') +
  geom_point(data = ret.h0.diff.pa08.c08.h10, aes(x = re.max, y = re.max.value, group = factor(epsilon)), size = 2, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  ylim(re.min, re.max * 1.3) +  #
  guides(colour = guide_legend(title = expression(epsilon), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev


ret.h0.pa05.c08.h10 <- subset(ret.h0, Pa == 0.5 & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa05.c08.h10 <- subset(ret.h0.diff, Pa == 0.5 & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa05.c08.h10 <- subset(ret.h0.diff.pa05.c08.h10, round(f, digits = 5) %in% round(seq(0.1,1.5,0.1),digits=5))

# Fig. 1
re.min <- min(ret.h0.pa05.c08.h10$re)
re.max <- max(ret.h0.pa05.c08.h10$re)
p.lev <- ggplot(data = ret.h0.pa05.c08.h10, mapping = aes(x = Pm, y = lev, colour = factor(f))) +  # , group = factor(f)
  geom_line(aes(x = Pm, y = re, colour = factor(f)), linetype = 3) +
  geom_line() +
  geom_point(data = ret.h0.diff.pa05.c08.h10, aes(x = re.equal.rs, y = re.equal.rs.value, group = factor(f)), size = 3, shape = 23, colour = 'black') +
  geom_point(data = ret.h0.diff.pa05.c08.h10, aes(x = re.max, y = re.max.value, group = factor(f)), size = 3, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  ylim(re.min, 15) +  # re.max * 1.5
  guides(colour = guide_legend(title = expression(f), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev
#grid_arrange_shared_legend(p.v, p.re, p.vc, p.rs, p.rho, p.lev)


ret.h0.pa0.c08.h10 <- subset(ret.h0, Pa == 0. & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa0.c08.h10 <- subset(ret.h0.diff, Pa == 0. & round(C, digits = 5) == 0.8 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa0.c08.h10 <- subset(ret.h0.diff.pa0.c08.h10, round(f, digits = 5) %in% round(seq(0.1,1.5,0.1),digits=5))

# Fig. 1
re.min <- min(ret.h0.pa0.c08.h10$re)
re.max <- max(ret.h0.pa0.c08.h10$re)
p.lev <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = lev, colour = factor(f))) +  # , group = factor(f)
  geom_line(aes(x = Pm, y = re, colour = factor(f)), linetype = 3) +
  geom_line() +
  geom_point(data = ret.h0.diff.pa0.c08.h10, aes(x = re.equal.rs, y = re.equal.rs.value, group = factor(f)), size = 3, shape = 23, colour = 'black') +
  geom_point(data = ret.h0.diff.pa0.c08.h10, aes(x = re.max, y = re.max.value, group = factor(f)), size = 3, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  ylim(re.min, 25) +  # re.max * 1.5
  guides(colour = guide_legend(title = expression(f), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev
#grid_arrange_shared_legend(p.v, p.re, p.vc, p.rs, p.rho, p.lev)

ret.h0.pa0.f05.h10 <- subset(ret.h0, Pa == 0. & round(f, digits = 5) == 0.5 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa0.f05.h10 <- subset(ret.h0.diff, Pa == 0. & round(f, digits = 5) == 0.5 & round(h, digits = 5) == 1.0)
ret.h0.diff.pa0.f05.h10 <- subset(ret.h0.diff.pa0.f05.h10, round(C, digits = 5) %in% round(seq(0.1,1.0,0.1),digits=5))
# Fig. 1
re.min <- min(ret.h0.pa0.f05.h10$re)
re.max <- max(ret.h0.pa0.f05.h10$re)
p.lev <- ggplot(data = ret.h0.pa0.f05.h10, mapping = aes(x = Pm, y = lev, colour = factor(C))) +  # , group = factor(f)
  geom_line(aes(x = Pm, y = re, colour = factor(C)), linetype = 3) +
  geom_line() +
  geom_point(data = ret.h0.diff.pa0.f05.h10, aes(x = re.equal.rs, y = re.equal.rs.value, group = factor(C)), size = 3, shape = 23, colour = 'black') +
  geom_point(data = ret.h0.diff.pa0.f05.h10, aes(x = re.max, y = re.max.value, group = factor(C)), size = 3, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  ylim(re.min, re.max*1.5) +  # re.max * 1.5
  guides(colour = guide_legend(title = expression(C), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev

ret.h0.c08.f05.h10 <- subset(ret.h0, round(C, digits = 5) == 0.8 & round(f, digits = 5) == 0.5 & round(h, digits = 5) == 1.0)
ret.h0.c08.f05.h10 <- subset(ret.h0.c08.f05.h10, round(Pa, digits = 5) %in% round(seq(0.,0.9,0.1),digits=5))
ret.h0.diff.c08.f05.h10 <- subset(ret.h0.diff, round(C, digits = 5) == 0.8 & round(f, digits = 5) == 0.5 & round(h, digits = 5) == 1.0)
ret.h0.diff.c08.f05.h10 <- subset(ret.h0.diff.c08.f05.h10, round(Pa, digits = 5) %in% round(seq(0.,0.9,0.1),digits=5))
# Fig. 1
re.min <- min(ret.h0.c08.f05.h10$re)
re.max <- max(ret.h0.c08.f05.h10$re)
p.lev <- ggplot(data = ret.h0.c08.f05.h10, mapping = aes(x = Pm, y = lev, colour = factor(Pa))) +  # , group = factor(f)
  geom_line(aes(x = Pm, y = re, colour = factor(Pa)), linetype = 3) +
  geom_line() +
  geom_point(data = ret.h0.diff.c08.f05.h10, aes(x = re.equal.rs, y = re.equal.rs.value, group = factor(Pa)), size = 3, shape = 23, colour = 'black') +
  geom_point(data = ret.h0.diff.c08.f05.h10, aes(x = re.max, y = re.max.value, group = factor(Pa)), size = 3, shape = 19, colour = 'black') +
  scale_colour_hue(h=c(90, 360)) +
  labs(x = expression(P[m]), y = expression(R)) +
  ylim(re.min, re.max*1.5) +  # re.max * 1.5
  guides(colour = guide_legend(title = expression(P[a]), title.theme = element_text(size=15,angle = 0),reverse = TRUE)) +
  theme_bw()
p.lev


require(ggplot2)
# p.v <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = v, colour = factor(f))) +  # , group = factor(f)
#   geom_line() +
#   scale_colour_hue("f", h=c(90, 360)) +
#   labs(x = expression(P[m]), y = expression(V)) + #, title = "(1)"
#   theme_bw()
# p.vc <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = vc, colour = factor(f))) +  # , group = factor(f)
#   geom_line() +
#   scale_colour_hue("f", h=c(90, 360)) +
#   labs(x = expression(P[m]), y = expression(V[c])) +
#   theme_bw()
# p.rho <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = rho, colour = factor(f))) +  # , group = factor(f)
#   geom_line() +
#   scale_colour_hue("f", h=c(90, 360)) +
#   labs(x = expression(P[m]), y = expression(rho)) +
#   theme_bw()
# p.re <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = re, colour = factor(f))) +  # , group = factor(f)
#   geom_line() +
#   geom_point(data = ret.h0.diff.pa0.c08.h10, aes(x = re.max, y = re.max.value, group = factor(f)), size = 2, shape = 19, colour = 'black') +
#   scale_colour_hue("f", h=c(90, 360)) +
#   labs(x = expression(P[m]), y = expression(Re)) +
#   theme_bw()
# p.rs <- ggplot(data = ret.h0.pa0.c08.h10, mapping = aes(x = Pm, y = re1deriv, colour = factor(f))) +  # , group = factor(f)
#   geom_line() +
#   scale_colour_hue("f", h=c(90, 360)) +
#   labs(x = expression(P[m]), y = expression(Rs)) +
#   theme_bw()


# Fig. 2
require(reshape2)
ret.h0.diff.pa0.h10 <- subset(ret.h0.diff, Pa == 0 & round(h, digits = 5) == 1.0)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="diff")
write.table(tmp, file = 'df2diffFG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="deep")
write.table(tmp, file = 'df2deepFG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="deep1")
write.table(tmp, file = 'df2deep1FG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="deepavg")
write.table(tmp, file = 'df2deepavgFG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="re.max")
write.table(tmp, file = 'df2remax.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa0.h10, C~epsilon, value.var="re.equal.rs")
write.table(tmp, file = 'df2reequalrs.csv', sep = ',', row.names = F, col.names = F)

ret.h0.diff.pa05.h10 <- subset(ret.h0.diff, Pa == 0.5 & round(h, digits = 5) == 1.0)
tmp <- acast(ret.h0.diff.pa05.h10, C~epsilon, value.var="diff")
write.table(tmp, file = 'df2diffpa05FG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa05.h10, C~epsilon, value.var="deep")
write.table(tmp, file = 'df2deeppa05FG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa05.h10, C~epsilon, value.var="deep1")
write.table(tmp, file = 'df2deep1pa05FG.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa05.h10, C~epsilon, value.var="deepavg")
#tmp[tmp==Inf] <- 0
write.table(tmp, file = 'df2deepavgpa05FG.csv', sep = ',', row.names = F, col.names = F)

ret.h0.diff.pa08.h10 <- subset(ret.h0.diff, Pa == 0.8 & round(h, digits = 5) == 1.0)
tmp <- acast(ret.h0.diff.pa08.h10, C~f, value.var="diff")
write.table(tmp, file = 'df2diffpa08.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa08.h10, C~f, value.var="deep")
write.table(tmp, file = 'df2deeppa08.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa08.h10, C~f, value.var="deep1")
write.table(tmp, file = 'df2deep1pa08.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa08.h10, C~f, value.var="deepavg")
write.table(tmp, file = 'df2deepavgpa08.csv', sep = ',', row.names = F, col.names = F)

ret.h0.diff.pa09.h10 <- subset(ret.h0.diff, Pa == 0.9 & round(h, digits = 5) == 1.0)
tmp <- acast(ret.h0.diff.pa09.h10, C~f, value.var="diff")
write.table(tmp, file = 'df2diffpa09.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa09.h10, C~f, value.var="deep")
write.table(tmp, file = 'df2deeppa09.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa09.h10, C~f, value.var="deep1")
write.table(tmp, file = 'df2deep1pa09.csv', sep = ',', row.names = F, col.names = F)
tmp <- acast(ret.h0.diff.pa09.h10, C~f, value.var="deepavg")
write.table(tmp, file = 'df2deepavgpa09.csv', sep = ',', row.names = F, col.names = F)

#------>  import 'df2**.csv' into Matlab, then use Matlab to figure.


# Fig. 3
tmp <- subset(ret.h0.v5, C == 0.5 & round(Pc, 5) != 0)
write.table(tmp[,c('Pa','Pm','Pc','lev')], file = 'tmp.csv', sep = ',',row.names = F, col.names = F)
c05 <- subset(ret.h0, C == 0.1 & round(Pc, 5) != 0)
write.table(c05[,c('Pa','Pm','Pc','lev')], file = 'c05.csv', sep = ',',row.names = F, col.names = F)
c05eps01 <- subset(ret.h0, C == 0.5 & round(epsilon, digits = 5) == 0.1 & h == 1.)
write.table(c05eps01[,c('Pa','Pm','Pc','lev')], file = 'c05eps01.csv', sep = ',',row.names = F, col.names = F)

require(ggplot2)
ggplot(df, aes(x = Pm, y = v.caculated.taylor3.at05, group = factor(Pa))) +
  geom_line(aes(colour=factor(Pa)), size = 1.) +
  geom_point(aes(colour=factor(Pa)), size = 1.) +
  labs(x = expression(P[m]), y = 'R', colour = expression(P[a])) +
  theme_bw()

plot(df$Pm, df$v.caculated - df$v.caculated.taylor3.at05, log = 'y')



plot.ggtern.h0 <- function(df, colname) {
  require(ggtern)
  # sequence: Pa --> Pc --> Pm
  ggtern(df,aes(Pa, Pc, Pm)) +
    geom_point(aes_string(fill=colname),color="black",shape=21,size=3, alpha = 0.5) +
    scale_fill_gradient(low="green",high="red") +
    theme(legend.position=c(0,1),legend.justification=c(0,1)) +
    labs(fill=colname, x = expression(P[a]), y = expression(P[c]), z = expression(P[m])) +
    theme_bw()
}
plot.ggtern.h0(df, 'lev')
tmp <- subset(df, Pa == 0)
plot(Pm, rho)
plot(Pm, re)
plot(Pm, rs)
plot(Pm, lev)
plot(Pm, re.sub.rs)
plot(Pm, re1deriv)



###############################################################################
#### Test correlation among stability measurements
###############################################################################
# test several measurements of stability
rho <- seq(from = 0, to = 60, by = 0.1)
stability <- sapply(rho, function(rho) {
  A = matrix(c(-1, -(1+rho)^2, 1, -sqrt(1)), nrow = 2)
  stability = get_stability(A)
})

stability <- as.data.frame(t(stability))
