library(bnlearn)
library(clue)
library(transport)
library(cluster)
library(seewave)
library(dplyr)

clustering = function(datos.ini, min_divisiones, metodo, max_clusters=15, distrib) {
  
  # Se transforma la matriz de datos en factores.
  datos.ini = as.data.frame(datos.ini)
  datosfac = as.data.frame(lapply(datos.ini, factor))
  
  # Se añade la variable de clase y se genera el primer modelo al azar.
  datos = cbind(datosfac, as.factor(rbinom(length(datosfac[, 1]), 1, 0.5)))
  names(datos)[length(names(datos))] = "C"
  
  # Se crea el grafo con la estructura del modelo Naive-Bayes.
  variables = colnames(datos)
  g = empty.graph(variables)
  for (i in 1:(length(variables) - 1)) {
    g = set.arc(g, "C", variables[i])
  }
  
  # Inicialización de variables
  c.valores = c(0, 1)
  i_ciclo_exterior = 1
  mejor.modelo = list()
  bic.mej.ant = -Inf
  bic.mej.post = 0
  bic.mejores = vector()
  
  # Inicio del bucle exterior
  while ((bic.mej.ant < bic.mej.post || i_ciclo_exterior <= min_divisiones) && 
         i_ciclo_exterior <= max_clusters) {
    
    # Inicialización de variables
    i_ciclo_interior = 1
    lista.modelos = list()
    bic.antiguo = -Inf
    bic.nuevo = score(g, data = datos, type = "bic")
    bic.grupos = bic.nuevo
    
    # Inicio del bucle interior
    while ((bic.antiguo < bic.nuevo || i_ciclo_interior <= 20) && i_ciclo_interior <= 150) {
      
      lista.modelos[[i_ciclo_interior]] = datos["C"]
      
      # Cálculo de probabilidades en el modelo
      bn = bn.fit(g, data = datos, method = "bayes", iss = 10)
      # P(Ci)
      pc = bn[[length(variables)]]$prob  
      # P(Xi = xi | Ci)
      px.c = matrix(nrow = length(datos[, 1]), ncol = length(variables) - 1) 
      # P(ci | x1, ..., xn)
      pc.x = matrix(nrow = length(datos[, 1]), ncol = length(levels(datos[, length(variables)])))  
      
      # Recorremos los resultados de bn y asignamos las probabilidades a las variables
      for (k in 1:length(c.valores)) {
        for (j in 1:(length(variables) - 1)) {
          px.c[, j] = bn[[j]]$prob[datos[, j], k]
        }
        pc.x[, k] = pc[k] * apply(px.c, 1, prod)
      }
      pc.x = pc.x / apply(pc.x, 1, sum)
      
      # Generación del modelo nuevo a partir de pc.x
      # usando pc.x como parametros de la multinomial
      rand = rmultinom(1, 1, pc.x[1, ])
      for (i in 2:length(datos[, 1])) {
        rand = cbind(rand, rmultinom(1, 1, pc.x[i, ]))
      }
      for (l in 1:length(c.valores)) {
        datos[which(rand[l, ] == 1), length(variables)] = l - 1
      }
      
      # Cálculo del BIC
      bic.antiguo = bic.nuevo
      bic.nuevo = score(g, data = datos, type = "bic")
      bic.grupos = append(bic.grupos, bic.nuevo)
      
      i_ciclo_interior = i_ciclo_interior + 1
    }
    # Final del bucle interior
    
    # Selección del mejor modelo
    bic.mejores[i_ciclo_exterior] = max(bic.grupos)
    mejor.modelo[[i_ciclo_exterior]] = lista.modelos[[min(which(bic.grupos == bic.mejores[i_ciclo_exterior]))]]
    
    dividir = -1
    # PARTE DONDE SE ELIJE EL PROXIMO CLUSTER A DIVIDIR
    if (metodo == "Wasserstein") {
      if(distrib == 'unif'){
      # EL PROXIMO GRUPO A DIVIDIR ES AQUEL CON DISTANCIA DE WASSERSTEIN
      # MAS CERCANA A LA DISTRIBUCION UNIFORME
        disim = vector(length = length(variables) - 1) # distancias de cada variable
        # suma de las distancias de las variables de cada grupo
        wasserstein.total = vector(length = length(c.valores)) 
        
        for (i in 1:length(c.valores)) {
          grupo = datos[which(datos$C == i - 1), ]
          for (q in 1:(length(variables) - 1)) {
            distrib_observada <- table(grupo[[variables[q]]]) / nrow(grupo)
            distrib_uniforme <- rep(1 / length(distrib_observada), length(distrib_observada))  # Uniforme
            
            disim[q] <- wasserstein1d(distrib_observada, distrib_uniforme)
          }
          wasserstein.total[i] <- sum(disim)
        }
        dividir <- which.min(wasserstein.total) - 1
        
      } else if(distrib == 'global'){
      # EL PROXIMO GRUPO A DIVIDIR ES AQUEL CON DISTANCIA DE WASSERSTEIN
      # MAS CERCANA A LA DISTRIBUCION GLOBAL DE CADA VARIABLE 
        disim = vector(length = length(variables) - 1) # distancias de cada variable
        # suma de las distancias de las variables de cada grupo
        wasserstein.total = vector(length = length(c.valores)) 
        
        # calcular la distribución global para cada variable
        distrib_global <- lapply(variables, function(var) {
          table(datos[[var]]) / nrow(datos)
        })
        
        for (i in 1:length(c.valores)) {
          grupo <- datos[which(datos$C == i - 1), ]
          for (q in 1:(length(variables) - 1)) {
            distrib_observada <- table(grupo[[variables[q]]]) / nrow(grupo)
            distrib_referencia <- distrib_global[[q]]
            
            disim[q] <- wasserstein1d(distrib_observada, distrib_referencia)
          }
          wasserstein.total[i] <- sum(disim)
        }
        dividir <- which.min(wasserstein.total) - 1
      }
    } else if (metodo == "KL") {
      if(distrib == 'unif'){
      # EL PROXIMO GRUPO A DIVIDIR ES AQUEL CON DISIMILITUD DE KL
      # MENOR A LA DISTRIBUCION UNIFORME 
        disim = vector(length = length(variables) - 1)
        kl.total = vector(length = length(c.valores))
        for (i in 1:length(c.valores)) {
          grupo = datos[which(datos$C == i - 1), ]
          for (q in 1:(length(variables) - 1)) {
            probs = table(grupo[, q]) / length(grupo[, q])
            disim[q] = kl.dist(probs, rep(1 / length(probs), length(probs)), base = exp(1))[[3]]
          }
          kl.total[i] = sum(disim)
        }
        dividir = which.min(kl.total) - 1
      } else if(metodo == 'global'){
      # EL PROXIMO GRUPO A DIVIDIR ES AQUEL CON DISIMILITUD DE KL
      # MENOR A LA DISTRIBUCION GLOBAL DE CADA VARIABLE 
        disim = vector(length = length(variables) - 1)
        kl.total = vector(length = length(c.valores))
        # calcular la distribución global para cada variable
        distrib_global <- lapply(variables, function(var) {
          table(datos[[var]]) / nrow(datos)
        })
        for (i in 1:length(c.valores)) {
          grupo = datos[which(datos$C == i - 1), ]
          for (q in 1:(length(variables) - 1)) {
            probs = table(grupo[, q]) / length(grupo[, q])
            disim[q] = kl.dist(probs, , base = exp(1))[[3]]
          }
          kl.total[i] = sum(disim)
        }
        dividir = which.min(kl.total) - 1
      }
    } else {
      stop("Por favor inserte un método válido")
    }
    
    # Adición del nuevo cluster
    c.valores = append(c.valores, length(c.valores))
    levels(datos$C) = c.valores
    
    datos[datos$C == dividir, "C"] = as.factor(sample(c(dividir, length(c.valores) - 1), length(datos[datos$C == dividir, "C"]), replace = TRUE))
    
    # Actualización de BIC de los mejores modelos
    if (i_ciclo_exterior > 1) {
      bic.mej.ant = bic.mejores[i_ciclo_exterior - 1]
      bic.mej.post = bic.mejores[i_ciclo_exterior]
    }
    
    i_ciclo_exterior = i_ciclo_exterior + 1
  }
  # Final del bucle exterior
  
  # Creación de la salida
  if(min_divisiones == max_clusters){
    mejor.mejor = mejor.modelo[[min_divisiones]]
  } else {
    mejor.mejor = mejor.modelo[[which(bic.mejores == max(bic.mejores))]]
  }
  salida = list(table(mejor.mejor), mejor.mejor, bic.mejores)
  return(salida)
}

