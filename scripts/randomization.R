##different methods for randomizing with correlation

local.coin.jitter <- function(vector, window){
  n <- length(vector)
  idx <- 1:n

  #because sample when idx[start:end] is length 1 defaults to constructing a vector of integers from 1:idx[start:end] when sampling, we have to go around this

  if (window==0){
    return(vector)
  }

  for (k in idx){
      start <- max(1, k-window)
      end <- min(k+window, n)
      j <- sample(idx[start:end], 1)
      swap <- vector[k]
      swip <- vector[j]
      vector[j] <- swap
      vector[k] <- swip
  }
  return(vector)
}

in.the.looking.glass <- function(vector, data, true, error){

  original <- vector
  n <- length(vector)
  idx <- 1:n

  true <- abs(true)

  if (true==1){
    return(vector)
  }

  cor <- cor(data$drugA, vector, method="spearman")

  loss <- abs(cor - true)

  window <- ceiling(n*0.01)
  increment <-  ceiling(n*0.1)

  counter <- 1

  if (cor < true){
    return(vector)
  }


  while (loss > error){

    vector <- local.coin.jitter(vector, window)
    cor <- cor(data$drugA, vector, method="spearman")
    loss <- abs(cor - true)

    if (counter > 9){
      counter <- 0
      window <- min(window+increment, n)
      vector <- original
    }


    if ((window >= n)|(window <= 0)){
      out <- find.correlated.vector(data, rho=true, param=40, stall=20, type="jitter",error=error)

      print("executed failsafe for simulation estimation")

      return(out)
    }

    counter <- counter + 1
  }

  return(vector)
}

local.coin.jitter2 <- function(vector, window){
  n <- length(vector)
  idx <- 1:n

  #because sample when idx[start:end] is length 1 defaults to constructing a vector of integers from 1:idx[start:end] when sampling, we have to go around this

  if (window==0){
    return(vector)
  }

  for (k in idx){
      start <- max(1, k-window)
      end <- min(k+window, n)
      j <- sample(idx[start:end], 1)
      swap <- vector[k]
      swip <- vector[j]
      vector[j] <- swap
      vector[k] <- swip
  }

  for (k in rev(idx)){
      start <- min(k+window, n)
      end <- max(1, k-window)
      j <- sample(idx[end:start], 1)
      swap <- vector[k]
      swip <- vector[j]
      vector[j] <- swap
      vector[k] <- swip
  }

  return(vector)
}

jitter.randomizer <- function(vector, window){
  idx <- 0
  store <- vector()
  while (length(vector) > window){
    start <- 1
    end <-  window
    temp <- 1:length(vector)
    subset <- temp[start:end]
    pick <- sample(subset, 1)
    idx <- idx + 1
    store[idx] <- vector[pick]
    vector <- vector[-pick]
  }
  if (window > 0){
  store[(idx+1):(idx+window)] <- sample(vector, length(vector), replace=FALSE)
  }
  return(store)
}

weighted.coin.randomizer <- function(vector, p){
  idx <- 1:length(vector)
  for (k in idx){
    out <- sample(c("stay", "switch"), 1, prob=c(1-p, p))
    if (out=="switch"){
      j <- sample(idx[-k], 1)
      swap <- vector[k]
      swip <- vector[j]
      vector[j] <- swap
      vector[k] <- swip
    }
  }
  return(vector)
}

#methods for figuring out how to shuffle a vector given a certain target correlation

find.correlated.vector <- function(data, rho, param, stall, type, error){
  vector <- data$drugB
  comp.vector <- data$drugA

  rho <- abs(rho)

  if (type=="double_walk"){
    size <- dim(data)[1]/2
    space <- seq(0, dim(data)[1], size)
    values <- lapply(space, function(x) local.coin.jitter2(vector,x))
    cor.values <- unlist(lapply(values, function(x) cor(x, comp.vector, method="spearman")))
  }

  if (type=="window.swap"){
    size <- dim(data)[1]/2
    space <- seq(0, dim(data)[1], size)
    values <- lapply(space, function(x) local.coin.jitter(vector,x))
    cor.values <- unlist(lapply(values, function(x) cor(x, comp.vector, method="spearman")))
  }

  if (type=="jitter"){
    size <- dim(data)[1]/2
    space <- seq(0, dim(data)[1], size)
    values <- lapply(space, function(x) jitter.randomizer(vector,x))
    cor.values <- unlist(lapply(values, function(x) cor(x, comp.vector, method="spearman")))
  }

  if (type=="coin"){
    space <- c(0,0.5,1)
    values <- lapply(space, function(x) weighted.coin.randomizer(vector, x))
    cor.values <- unlist(lapply(values, function(x) cor(x, comp.vector, method="spearman")))
  }

  if (rho < min(cor.values) | rho > max(cor.values)){
    init.cor <- cor(vector, comp.vector, method="spearman")
    exit.condition <- rep(1,30)
    idx.count <- 1

    if (rho > max(cor.values)){
      coin.p <- 0.0000001
      jitter.exit <- TRUE
    } else{
      coin.p <- 1
      jitter.exit <- FALSE
    }

    if (type=="coin"){
      while (abs(init.cor - rho) > error){
        vector <- weighted.coin.randomizer(vector, coin.p)
        init.cor <- cor(vector, comp.vector, method="spearman")

        exit.condition[idx.count] <- abs(init.cor - rho)
        idx.count <- idx.count + 1
        if ((max(exit.condition) - min(exit.condition)) < 0.01){
            return(vector)

            #delete this line after finding weights
            #return(list(vector=vector,weight=1))
        }
      }
    }

    if ((type=="jitter")|(type=="window.swap")|(type=="double_walk")){
      if (jitter.exit){
          return(vector)
      }

      while (abs(init.cor - rho) > error){
        vector <- jitter.randomizer(vector,dim(data)[1])
        init.cor <- cor(vector, comp.vector, method="spearman")

        exit.condition[idx.count] <- abs(init.cor - rho)
        idx.count <- idx.count + 1

        if ((max(exit.condition) - min(exit.condition)) < 0.01){
            return(vector)
        }
      }
    }

    return(vector)
    #delete this line after finding weights
    #return(list(vector=vector,weight=1))
  }

  loss <- min(abs(rho-cor.values))
  randomized.vector <- values[[which.min(abs(rho-cor.values))]]

  if (loss < error){
    #delete after get weights
    #return(list(vector=randomized.vector,weight=0))

    return(randomized.vector)
  }

  vector <- iterative.narrower(loss=loss, error=error, space=space, rho=rho, cor.values=cor.values, type=type, margin=param, input=vector,stall=stall, comp.vector=comp.vector)

  return(vector)

}

iterative.narrower <- function(loss, error, space, rho, cor.values, type, margin, input, stall, comp.vector){

  cache <- rep(1, length=stall)
  counter <- 1

  while(loss > error){
    interval <- space[order(abs(rho-cor.values))][1:2]

    if ((type=="jitter")|(type=="window.swap")|(type=="double_walk")){
      size <- ceiling((max(interval) - min(interval))/2)
      start <- min(interval)-margin
      end <- max(interval)+margin
      if (start < 0){
        start <- 0
      }
      if (end > length(comp.vector)){
        end <- length(comp.vector)
      }
    }


    if (type=="coin"){
      size <- (max(interval) - min(interval))/2
      start <- min(interval) - margin
      end <- max(interval) + margin
      if (start < 0){
        start <- 0
      }
      if (end > 1){
        end <- 1
      }
    }

    space <- seq(start, end, size)

    if (type=="jitter"){
      values <- lapply(space, function(x) jitter.randomizer(input, x))
    }

    if (type=="window.swap"){
      values <- lapply(space, function(x) local.coin.jitter(input, x))
    }

    if (type=="double_walk"){
      values <- lapply(space, function(x) local.coin.jitter(input, x))
    }

    if (type=="coin"){
      values <- lapply(space, function(x) weighted.coin.randomizer(input, x))

    }

    cor.values <- unlist(lapply(values, function(x) cor(x, comp.vector, method="spearman")))
    loss <- min(abs(rho - cor.values))
    output <- values[[which.min(abs(rho - cor.values))]]

    #delete this line after getting weight
    #temp.solution <- space[which.min(abs(rho - cor.values))]

    if ((counter %% stall)==0){
      cache[stall] <- loss
    } else{
      cache[(counter %% stall)] <- loss
    }

    counter <- counter + 1

    if ((max(cache) - min(cache)) < 0.08){
      return(output)

      #delete after weight retreived
      #return(list(vector=output, weight=temp.solution))
    }

  }
  return(output)
  #delete after weight retreived
  #return(list(vector=output, weight=temp.solution))
}
