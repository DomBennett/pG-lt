## No copyright, no warranty
## Dominic John Bennett & Junying Lim
## Last update: 21/07/2013

## Libraries
library(picante) # Phylogenetic analytical tools
library(caper)
library(lme4)
library(stringr) #String manipulation tools
library(RCurl) ##GET and POST to URL queries functionality
library(apTreeshape) # randI

## Dom's Functions:
deg2dec <- function(data) {
  # Takes 'coordinate strings' in degrees (D M (S)) and converts to decimals.
  # Expects coordinate string such: ##'##'##'
  #   - where: # denotes any number of any length, ' denotes any non-numeric(s).
  # Function will stop if more than 3 numbers in coordinate string
  #
  # Args:
  #   data: a matrix or dataframe of coordinates (Latitude and Longitude).
  #
  # Returns:
  #   A dataframe of converted coordinates.
  
  # input checking
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("Input must be a matrix or dataframe, such:\n
         Latitude     Longitude\n
         51 28\'46.92\"  0 0\'29.47\"")
  }
  if(ncol(data)!= 2) {
    stop("Input must have two columns, such:\n
         Latitude     Longitude\n
         51 28\'46.92\"  0 0\'29.47\"")
  }
  two.numbers <- apply(data, 1:2, function(x) grepl("([^0-9][0-9]|[0-9][^0-9])", x))
  if(!all(two.numbers)) {
    print("Invalid data entries:")
    print(data[!two.numbers])
    stop("Input must be coordinate character strings with at least 2 numbers, e.g.:\n
         Latitude     Longitude\n
         51 28\'46.92\"  0 0\'29.47\"")
  }
  
  # vector inputs
  x <- data[ ,1] #latitude is x
  y <- data[ ,2] #longitude is y
  
  # find start and end of numbers in coordinate strings
  index.x1 <- gregexpr("[0-9\\-\\.]+", x)
  index.x2 <- gregexpr("([0-9\\.][^0-9\\.])|[0-9]$", x)
  index.y1 <- gregexpr("[0-9\\-\\.]+", y)
  index.y2 <- gregexpr("([0-9\\.][^0-9\\.])|[0-9]$", y)
  
  # vector outputs
  decimal.x <- rep(NA, nrow(data))
  decimal.y <- rep(NA, nrow(data))
  
  # loop through each row of coordinate string:
  #   - take the first number in string, add the second number/60
  #     and then the third/3600
  #    - stop if more than 3 numbers found.
  error.message1 <- "\nMore than 3 numbers in:\n data["
  error.message2 <- "\nTry checking for non-numeric characters in numbers e.g.:\n
  Here, there is an extra space between the 9 and 2 -- 51 28\'46.9 2\""
  for(i in 1:nrow(data)) {
    decimal.x[i] <- as.numeric(substr(x[i], index.x1[[i]][1], index.x2[[i]][1]))
    for(j in 2:length(index.x1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 1, ",", i,"] -- ",x[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.x <- as.numeric(substr(x[i], index.x1[[i]][j], index.x2[[i]][j]))
      decimal.x[i] <- decimal.x[i] + (next.x/division)
    }
    decimal.y[i] <- as.numeric(substr(y[i], index.y1[[i]][1], index.y2[[i]][1]))
    for(j in 2:length(index.y1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 2, ",", i,"] -- ",y[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.y <- as.numeric(substr(y[i], index.y1[[i]][j], index.y2[[i]][j]))
      decimal.y[i] <- decimal.y[i] + (next.y/division)
    }
  }
  
  # control for S and W
  decimal.x <- decimal.x * ifelse(grepl("(S|s)", x), -1, 1)
  decimal.y <- decimal.y * ifelse(grepl("(W|w)", y), -1, 1)
  
  # generate output
  decimals <- data.frame(Latitude = decimal.x, Longitude = decimal.y)
  return(decimals)
}

genNullDist <- function(phylo, comm.data, htypes, null = 0, nrands = 2000,
                        metric = "pse", abundance.weighted = TRUE, pd.type = 1,
                        pd.min = 0) {
  # Generate null distributions of differences for the specified phylogenetic
  #  metric between habitat type 1 and habitat type 2 (1 - 2) for community data
  #  using null models as described by Helmus et al 2010.
  #
  # Args:
  #  phylo: phylo object of the community data
  #  comm.data: community matrix (species columns, sites rows)
  #  htypes: vector of habitat types for each site (only 2 types allowed!)
  #  null: number of null distribution. Null = 0, randomisation test, Nulls 1 to 5
  #   null distributions as described by Helmus et al 2010
  #  nrands: number of iterations
  #  metric: either pse, psv, mpd or pd at this stage
  #
  # Return:
  #  vector of null differences
  #
  # TODO(05/07/2013): This code is HORRENDOUSLY long! I've done this to avoid if
  #  statements in for loops, but there must be a better way of doing it!
  randPrev <- function(mat) {
    ## Randomise prevalences
    spp.names <- colnames(mat)
    out <- t(apply(mat, 1, sample))
    colnames(out) <- spp.names
    return (out)
  }
  randRich <- function(mat) {
    ## Randomise richnesses
    spp.names <- colnames(mat)
    out <- apply(mat, 2, sample)
    colnames(out) <- spp.names
    return (out)
  }
  randomise <- function(mat) {
    ## Randomise both
    mat <- randPrev(mat)
    mat <- randRich(mat)
    return (mat)
  }
  if (metric == "psv") {
    null0 <- function(null.dist) {
      ## Randomise richness and prevalence
      for (i in 1:nrands) {
        temp.data <- randomise(comm.data)
        temp.phymet <- psv(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null1 <- function(null.dist) {
      ## Randomise site richnesses
      for (i in 1:nrands) {
        temp.data <- randRich(comm.data)
        temp.phymet <- psv(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null2 <- function(null.dist) {
      ## Randomise species prevalences
      for (i in 1:nrands) {
        temp.data <- randPrev(comm.data)
        temp.phymet <- psv(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null3 <- function(null.dist) {
      ## Randomise within the habitat types
      u.htypes <- unique(htypes)
      data.htype.1 <- comm.data[htypes == u.htypes[1], ]
      data.htype.2 <- comm.data[htypes == u.htypes[2], ]
      for (i in 1:nrands) {
        temp.data <- rbind(randRich(data.htype.1), randRich(data.htype.2))
        temp.phymet <- psv(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null45 <- function(null.dist, data.rand, data.frozen) {
      ## 4 -- Maintain the species that are gained
      ## 5 -- Maintain the species that are lost
      for (i in 1:nrands) {
        temp.data <- cbind(randPrev(data.rand), data.frozen)
        temp.phymet <- psv(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null.dist <- rep(NA, nrands)
    if (null == 0) {
      null.dist <- null0(null.dist)
    } else if (null == 1) {
      null.dist <- null1(null.dist)
    } else if (null == 2) {
      null.dist <- null2(null.dist)
    } else if (null == 3) {
      null.dist <- null3(null.dist)
    } else if (null == 4 | null == 5){
      u.htypes <- unique(htypes)
      diffs <- colSums(comm.data[htypes == u.htypes[1], ]) -
        colSums(comm.data[htypes == u.htypes[2], ])
      if (null == 4) {
        gains <- ifelse(diffs < 0, 1, 0)
        data.rand <- comm.data[ , gains == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , gains == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- psv(comm.data, phylo)[ ,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      } else {
        losses <- ifelse(diffs > 0, 1, 0)
        data.rand <- comm.data[ , losses == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , losses == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- psv(data.frozen, phylo)[ ,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      }
    } else {
      stop("Null must be a number 0-5")
    }
  } else if (metric == "pse") {
    null0 <- function(null.dist) {
      ## Randomise richness and prevalence
      for (i in 1:nrands) {
        temp.data <- randomise(comm.data)
        temp.phymet <- pse(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null1 <- function(null.dist) {
      ## Randomise site richnesses
      for (i in 1:nrands) {
        temp.data <- randRich(comm.data)
        temp.phymet <- pse(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null2 <- function(null.dist) {
      ## Randomise species prevalences
      for (i in 1:nrands) {
        temp.data <- randPrev(comm.data)
        temp.phymet <- pse(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null3 <- function(null.dist) {
      ## Randomise within the habitat types
      u.htypes <- unique(htypes)
      data.htype.1 <- comm.data[htypes == u.htypes[1], ]
      data.htype.2 <- comm.data[htypes == u.htypes[2], ]
      for (i in 1:nrands) {
        temp.data <- rbind(randRich(data.htype.1), randRich(data.htype.2))
        temp.phymet <- pse(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null45 <- function(null.dist, data.rand, data.frozen) {
      ## 4 -- Maintain the species that are gained
      ## 5 -- Maintain the species that are lost
      for (i in 1:nrands) {
        temp.data <- cbind(randPrev(data.rand), data.frozen)
        temp.phymet <- pse(temp.data, phylo)[ ,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null.dist <- rep(NA, nrands)
    if (null == 0) {
      null.dist <- null0(null.dist)
    } else if (null == 1) {
      null.dist <- null1(null.dist)
    } else if (null == 2) {
      null.dist <- null2(null.dist)
    } else if (null == 3) {
      null.dist <- null3(null.dist)
    } else if (null == 4 | null == 5){
      u.htypes <- unique(htypes)
      diffs <- colSums(comm.data[htypes == u.htypes[1], ]) -
        colSums(comm.data[htypes == u.htypes[2], ])
      if (null == 4) {
        gains <- ifelse(diffs < 0, 1, 0)
        data.rand <- comm.data[ , gains == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , gains == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- pse(comm.data, phylo)[ ,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      } else {
        losses <- ifelse(diffs > 0, 1, 0)
        data.rand <- comm.data[ , losses == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , losses == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- pse(data.frozen, phylo)[ ,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      }
    } else {
      stop("Null must be a number 0-5")
    }
  } else if (metric == "mpd"){
    null0 <- function(null.dist) {
      ## Randomise richness and prevalence
      for (i in 1:nrands) {
        temp.data <- randomise(comm.data)
        temp.phymet <- mpd(temp.data, cophenetic(phylo),
                           abundance.weighted = abundance.weighted)
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null1 <- function(null.dist) {
      ## Randomise site richnesses
      for (i in 1:nrands) {
        temp.data <- randRich(comm.data)
        temp.phymet <- mpd(temp.data, cophenetic(phylo),
                           abundance.weighted = abundance.weighted)
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null2 <- function(null.dist) {
      ## Randomise species prevalences
      for (i in 1:nrands) {
        temp.data <- randPrev(comm.data)
        temp.phymet <- mpd(temp.data, cophenetic(phylo),
                           abundance.weighted = abundance.weighted)
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null3 <- function(null.dist) {
      ## Randomise within the habitat types
      u.htypes <- unique(htypes)
      data.htype.1 <- comm.data[htypes == u.htypes[1], ]
      data.htype.2 <- comm.data[htypes == u.htypes[2], ]
      for (i in 1:nrands) {
        temp.data <- rbind(randRich(data.htype.1), randRich(data.htype.2))
        temp.phymet <- mpd(temp.data, cophenetic(phylo),
                           abundance.weighted = abundance.weighted)
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null45 <- function(null.dist, data.rand, data.frozen) {
      ## 4 -- Maintain the species that are gained
      ## 5 -- Maintain the species that are lost
      for (i in 1:nrands) {
        temp.data <- cbind(randPrev(data.rand), data.frozen)
        temp.phymet <- mpd(temp.data, cophenetic(phylo),
                           abundance.weighted = abundance.weighted)
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null.dist <- rep(NA, nrands)
    if (null == 0) {
      null.dist <- null0(null.dist)
    } else if (null == 1) {
      null.dist <- null1(null.dist)
    } else if (null == 2) {
      null.dist <- null2(null.dist)
    } else if (null == 3) {
      null.dist <- null3(null.dist)
    } else if (null == 4 | null == 5){
      u.htypes <- unique(htypes)
      diffs <- colSums(comm.data[htypes == u.htypes[1], ]) -
        colSums(comm.data[htypes == u.htypes[2], ])
      if (null == 4) {
        gains <- ifelse(diffs < 0, 1, 0)
        data.rand <- comm.data[ , gains == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , gains == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- mpd(temp.data, cophenetic(phylo),
                             abundance.weighted = abundance.weighted)
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      } else {
        losses <- ifelse(diffs > 0, 1, 0)
        data.rand <- comm.data[ , losses == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , losses == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- mpd(temp.data, cophenetic(phylo),
                             abundance.weighted = abundance.weighted)
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      }
    } else {
      stop("Null must be a number 0-5")
    }
  } else if (metric == "pd") {
    null0 <- function(null.dist) {
      ## Randomise richness and prevalence
      for (i in 1:nrands) {
        temp.data <- randomise(comm.data)
        temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null1 <- function(null.dist) {
      ## Randomise site richnesses
      for (i in 1:nrands) {
        temp.data <- randRich(comm.data)
        temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null2 <- function(null.dist) {
      ## Randomise species prevalences
      for (i in 1:nrands) {
        temp.data <- randPrev(comm.data)
        temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null3 <- function(null.dist) {
      ## Randomise within the habitat types
      u.htypes <- unique(htypes)
      data.htype.1 <- comm.data[htypes == u.htypes[1], ]
      data.htype.2 <- comm.data[htypes == u.htypes[2], ]
      for (i in 1:nrands) {
        temp.data <- rbind(randRich(data.htype.1), randRich(data.htype.2))
        temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null45 <- function(null.dist, data.rand, data.frozen) {
      ## 4 -- Maintain the species that are gained
      ## 5 -- Maintain the species that are lost
      for (i in 1:nrands) {
        temp.data <- cbind(randPrev(data.rand), data.frozen)
        temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
        habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
        null.dist[i] <- habitat.means[1] - habitat.means[2]
      }
      return (null.dist)
    }
    null.dist <- rep(NA, nrands)
    if (null == 0) {
      null.dist <- null0(null.dist)
    } else if (null == 1) {
      null.dist <- null1(null.dist)
    } else if (null == 2) {
      null.dist <- null2(null.dist)
    } else if (null == 3) {
      null.dist <- null3(null.dist)
    } else if (null == 4 | null == 5){
      u.htypes <- unique(htypes)
      diffs <- colSums(comm.data[htypes == u.htypes[1], ]) -
        colSums(comm.data[htypes == u.htypes[2], ])
      if (null == 4) {
        gains <- ifelse(diffs < 0, 1, 0)
        data.rand <- comm.data[ , gains == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , gains == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      } else {
        losses <- ifelse(diffs > 0, 1, 0)
        data.rand <- comm.data[ , losses == 0]
        spp.names <- colnames(data.rand)
        data.rand <- as.matrix(data.rand)
        colnames(data.rand) <- spp.names
        data.frozen <- comm.data[ , losses == 1]
        spp.names <- colnames(data.frozen)
        data.frozen <- as.matrix(data.frozen)
        colnames(data.frozen) <- spp.names
        if (length(data.rand) <= nrow(comm.data)) {
          # if 1 or fewer species decrease -- no randomisation can be performed
          temp.phymet <- commPD(phylo, temp.data, type = pd.type, min.spp = pd.min)[,1]
          habitat.means <- tapply(temp.phymet, htypes, mean, na.rm = TRUE)
          diff <- habitat.means[1] - habitat.means[2]
          null.dist <- rep(diff, nrands)
        } else if (length(data.frozen) <= nrow(comm.data)){
          # if all - 1 or more species increase -- no randomisation can be performed
          null.dist <- null2(null.dist)
        } else {
          null.dist <- null45(null.dist, data.rand, data.frozen)
        }
      }
    } else {
      stop("Null must be a number 0-5")
    }
  } else {
    stop("Invalid metric provided. Must be either psv, pse, pd or mpd.")
  }
  return (null.dist)
}

meanPhylo <- function(phylo.dist, phylo.size = 5, min.phylos = 2) {
  # Take a distribution of phylogenies and return the mean phylogeny, using mean
  #  branch distances and the neighbour-joining algorithm
  #
  # Args:
  #  phylo.dist: a list of phylo objects with shared tips
  #  phylo.size: the minimum number of tips for a mean to be calculated
  #  min.phylos: the minimum number of phylogenies for a mean to be calculated
  #
  # Returns:
  #  a phylo object
  phylo.sizes <- sapply(phylo.dist, function(x) length(x$tip.label))
  if (max(phylo.sizes) < phylo.size) stop("All phylogenies have fewer tips than phylo.size")
  # I want to have the phylogenies 'overlapping' so drop ones that are too big
  #  to avoid having species that are represented by few phylogenies
  max.size <- mean(phylo.sizes) + mean(phylo.sizes) * 0.75
  min.size <- mean(phylo.sizes) - mean(phylo.sizes) * 0.25
  keep <- phylo.sizes > min.size & phylo.sizes < max.size 
  if (sum(keep) < min.phylos) stop("Fewer than min.phylos phylogenies of appropriate size")
  phylo.dist <- phylo.dist[keep]
  phylo.mats <- list()
  species.vec <- vector()
  # generate distance matrices for each phylogeny
  for (i in 1:length(phylo.dist)) {
    phylo.mats <- c(phylo.mats, list(cophenetic.phylo(phylo.dist[[i]])))
    species.vec <- c(species.vec, phylo.dist[[i]]$tip.label)
  }
  species.vec <- unique(species.vec)
  combs <- combn(species.vec, 2)
  comb.dists <- vector("list", ncol(combs))
  # extract distances from matrices for each combinate
  for (i in 1:length(phylo.mats)) {
    for (j in 1:ncol(combs)) {
      c <- combs[1, j] == colnames(phylo.mats[[i]])
      r <- combs[2, j] == rownames(phylo.mats[[i]])
      if (any(c) & any(r)) {
        comb.dists[[j]] <- c(comb.dists[[j]], phylo.mats[[i]][c, r])
      }
    }
  }
  mean.mat <- matrix(rep(0, length(species.vec)^2),
                     nrow = length(species.vec), ncol = length(species.vec))
  rownames(mean.mat) <- colnames(mean.mat) <- species.vec
  # take means and put into a matrix
  for (i in 1:ncol(combs)) {
    temp.dist <- mean(comb.dists[[i]])
    sp1.pos <- species.vec == combs[1, i]
    sp2.pos <- species.vec == combs[2, i]
    mean.mat[sp1.pos, sp2.pos] <- temp.dist
    mean.mat[sp2.pos, sp1.pos] <- temp.dist
  }
  return (nj(mean.mat))
}

extractEdges <- function(phylo, taxa, type = 1) {
  # Extract edges from a phylo object using 1 of 3 methods
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  taxa: vector of taxon names
  #  type:
  #     1 -- phylogeny consisting solely of the taxa, default
  #     2 -- edges from taxon tips to terminal node
  #     3 -- edges unique to taxa
  #
  # Return:
  #  vector of edges
  # TODO(01/07/2013): this may be more achievable with a vegan matrix
  if (!type %in% c(1,2,3)) {
    stop("Type must be an integer: 1, 2 or 3.")
  }
  if (!is.vector(taxa) | !is.character(taxa)) {
    stop("Invalid or no taxa given.")
  }
  if (length(taxa) == length (phylo$tip.label)){
    return(phylo$edge)
  }
  if (type == 1 & length (taxa) == 1){
    stop("length(taxa) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  edges <- match (match (taxa, phylo$tip.label), phylo$edge[,2])
  end.nodes <- phylo$edge[edges, 1]
  term.node <- length (phylo$tip.label) + 1
  if (all(end.nodes %in% term.node)) {
    return(edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (phylo$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, phylo$edge[,2])
          end.node <- phylo$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, phylo$edge[,2])
        end.node <- phylo$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
      }
    }
    return (edges)
  }
}

commPD <- function(phylo, comm.data, type = 1, min.spp = 2,
                   taxon.names = colnames(comm.data)) {
  # Calculate Faith's PD for a community given a community phylogeny and a community
  #  matrix. Depends on extractEdges().
  #
  # Args:
  #  phylo: community phylogeny (ape class)
  #  comm.data: matrix of community data, species as cols sites as rows
  #  type: which branches to include for calculating PD, either 1, 2 or 3,
  #   default 2 -- see extractEdges()
  #  min.spp: minimum number of species to include for calculating PD, default 2
  #  taxon.names: if taxon.names are not the column names of comm.data specify
  #   here, else ignore.
  #
  # Return:
  # vector of PDs
  if (type == 1 & min.spp < 2) {
    stop("Cannot compute type 2 PD for fewer than 2 species")
  }
  calcPD <- function (row) {
    if (sum(row) > 1){
      taxa <- taxon.names[as.logical(row)]
      edges <- extractEdges(phylo, taxa, type)
      return (sum(phylo$edge.length[edges]))
    } else {
      return (0)
    }
  }
  # add site names if none
  if (is.null(rownames(comm.data))) {
    rownames(comm.data) <- 1:nrow(comm.data)
  }
  # convert to incidence
  presences <- comm.data[comm.data > 0]
  comm.data[comm.data > 0] <- rep(1, length(presences))
  # drop sites with too few species
  row.sums <- rowSums(comm.data) >= min.spp
  site.names <- rownames(comm.data)[row.sums]
  comm.data <- comm.data[row.sums, ]
  # calculate PD by site
  pds <- as.matrix(apply(comm.data, 1, calcPD))
  rownames(pds) <- site.names
  return (pds)
}

plotComm <- function(comm.data, phylo, groups = rep(1, nrow(comm.data)),
                     no.margin = TRUE, ...){
  # Plot community data on community phylogeny to visualise
  #  differences in community structure. Use colours to distinguish
  #  site groups and alpha to distinguish abundances (works like rainbow())
  #
  # Args:
  #  comm.data: community data matrix (cols taxa, rows sites)
  #  phylo: community phylogeny
  #  groups: site groups
  #
  # Return:
  #  a matrix of community data
  # ultrametricize tree FIRST
  phylo<- compute.brlen(phylo, method="Grafen")
  # plot phylogeny, allow space for points
  edges <- extractEdges(phylo, phylo$tip.label[1], type = 2)
  phylo$edge.length <- phylo$edge.length/sum(phylo$edge.length[edges])
  # make all phylos the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(comm.data) < 20) {
    variable.max <- 1 + nrow(comm.data)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(comm.data)/10
    spacing.i <- 0.1 - 1/nrow(comm.data)
    # some stupid f***ing bizarre eqn! Why is it not linear!?!?!?!
    spacing.start <- 0.5 + spacing.i
  }
  plot(phylo, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max), ...)
  
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- comm.data/max(comm.data)
  
  # loop init
  ntips <- length(phylo$tip.label)
  spacing <- spacing.start
  group <- groups[1]
  j <- 1
  
  # loop through sites and plot points for present species
  for(i in 1:nrow(comm.data)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(comm.data[i,])
    taxa <- phylo$tip.label[pull]
    abunds <- alphas[i, pull]
    tiplabels(tip = match(taxa, phylo$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + spacing.i
    group <- groups[i]
  }
}

randCommData <- function(phylo, nsites, nspp, pa = TRUE, lam = 0.1){
  # Generate random community data for testing community phylogenetic
  #  methods
  #
  # Args:
  #   phylo: phylogeny
  #   nsites: number of sites
  #   nspp: number of species for each site
  #   pa: presence/absence data?
  #   lam: lambda for poisson distribution if pa is False
  #
  # Return:
  #   a matrix of community data
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  if (pa) {
    for (i in 1:nsites) {
      output[i, ] <- sample(c(rep(1, nspp), rep(0, ntips - nspp)))
    }
  } else {
    for (i in 1:nsites) {
      output[i, ] <- sample(rpois(ntips, lam))
    }
  }
  return(output)
}

split0 <- function(r = c(1, 10) , n = 2) {
  # evenly split a range into n
  #
  # args:
  #   r: range to be split, must be numerical, can be vector.
  #   n: times to split r
  #
  # returns a vector of splits
  output <- rep(NA, n)
  output[1] <- 0
  division <- 1/(n-1)
  for (i in 2:n) {
    output[i] <- output[i-1] + division
  }
  return(round((output * (max(r) - min(r))) + min(r)))
}

genCommData <- function(phylo, focal, fact = 1, mean.incid, mean.abun = FALSE,
                        nsites = 1) {
  # Generate clustered/overdispersed data for testing community
  #  phylogenetic analyses. Return either incidence or abundance data.
  #
  # Args:
  #  phylo - phylo.object, on which the community data will be based
  #  focal - numeric index, indicating which tip to perform cluster/dispersion
  #  fact - scaling coefficient: larger the number the more pronounced the
  #   effect by a power law (hence 0 not allowed)
  #  mean.incid - the mean incidence of species in the community
  #  mean.abun - the mean abundance per site, if given output will be abundances
  #  nsites - number of sites
  #
  # Return:
  #  a matrix of community data
  invertVector <- function(dists) {
    # for reversing the probs for overdispersion
    u <- sort(unique(dists), TRUE)
    s <- sort(u)
    probs <- rep(NA, length(dists))
    for (i in 1:length(u)) {
      probs[u[i] == dists] <- s[i]
    }
    return (probs)
  }
  genAbuns <- function(row) {
    # for generating abundances for each row
    out.row <- rep(0, ntips)
    temp.probs <- probs
    temp.probs[row < 1] <- 0
    abundance <- abs(ceiling(rnorm(1, mean = mean.abun)))
    if (abundance == 0) {
      return (out.row)
    } else {
      abuns <- sample(1:ntips, abundance, prob = temp.probs, replace = TRUE)
      abuns <- table(abuns)
      out.row[as.numeric(names(abuns))] <- abuns
      return (out.row)
    }
  }
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  pd.dists <- cophenetic.phylo(phylo)
  focal.dists <- pd.dists[ , focal] + 1 # avoid 0
  if (fact > 0) focal.dists <- invertVector(focal.dists)
  probs <- focal.dists^abs(fact)
  probs <- probs/sum(probs)
  for (i in 1:nsites) {
    incidence <- abs(ceiling(rnorm(1, mean = mean.incid))) # avoid negative numbers
    output[i, ] <- ifelse(1:ntips %in% sample(1:ntips, incidence,
                                              prob = probs), 1, 0)
  }
  if (mean.abun != FALSE) {
    for (i in 1:nsites) {
      output[i, ] <- genAbuns(output[i, ])
    }
  }
  return (output)
}

evenCommData <- function(phylo, nsites, nspp) {
  # Generate evenly distributed community phylogenetic data based on a
  #  phylogeny. Takes the phylogeny, calculates cummulative PD from random
  #  taxon, cuts at even intervals based on nspp to maximise distance.
  #
  # Args:
  #  phylo - phylo.object, on which community will be based
  #  nsites - number of sites
  #  nspp - number of species
  #
  # Depends:
  #  split0() - for evenly cutting up the cumsum.pd.dists
  #
  # Return:
  #  a matrix of community data 
  ntips <- length(phylo$tip.label)
  output <- matrix(rep(0, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- phylo$tip.label
  pd.dists <- cophenetic.phylo(phylo)
  for (i in 1:nsites) {
    focal.pd.dists <- pd.dists[ , sample(1:ntips, 1)]
    # select focal taxon at random
    splits <- split0(1:sum(focal.pd.dists), nspp)
    cumsum.pd.dists <- cumsum(focal.pd.dists)
    for (j in 1:nspp) {
      pull <- which(abs(cumsum.pd.dists - splits[j])
                    == min(abs(cumsum.pd.dists - splits[j])))
      output[i, pull] <- 1
    }
  }
  return(output)
}

phyloRarefy <- function(comm.data, phylo, samp, metric = 'PD', nrands = 2000,
                        type = 1) {
  # Monte-Carlo method for rarefiying community data using common phylogenetic
  #  metrics
  #
  # Args:
  #  comm.data: community matrix, species as cols and sites as rows
  #  phylo: phylogeny (ape class)
  #  samp: sample size (minimum incidence of abundance)
  #  metric: either PD, PSV, PSE, PSR, PSC or PSD
  #  nrands: the number of iterations
  #  types: if PD, type of PD to calculate
  #
  # Returns:
  #  matrix of rarefied PD by site and standard error
  if (ncol(comm.data) != length(phylo$tip.label)) {
    stop("Community data and phylogeny are different sizes.")
  } 
  if (!any(phylo$tip.label %in% colnames(comm.data))) {
    stop("Phylogeny tip names and community data names do not match.")
  }
  if (any(rowSums(comm.data) < samp)) {
    stop("Some sites have incidence/abundance less than samp")
  }
  # is it incidence data or abundance?
  if (any(comm.data > 1)) {
    abun <- TRUE
  } else {
    abun <- FALSE
  }
  # define internal functions for rapid vectorization
  # a lot of function definitions but avoids if statements for speed!!
  pluck <- function(row) {
    # this takes a community matrix and randomly distributes the occurences
    # of species using samp works for both incidence and abundance
    samples <- table(sample(1:length(row), samp, prob = row, replace = abun))
    sampled.row <- rep(0, length(row))
    sampled.row[as.numeric(names(samples))] <- samples
    return (sampled.row)
  } 
  PDrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(commPD(phylo, plucked, type, samp)))
  }
  PSVrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psv(plucked, phylo)[ ,1]))
  }
  PSErarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(pse(plucked, phylo)[ ,1]))
  }
  PSRrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psr(plucked, phylo)[ ,1]))
  }
  PSCrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psc(plucked, phylo)[ ,1]))
  }
  PSDrarefy <- function(x) {
    plucked <- t(apply(mat, 1, pluck))
    colnames(plucked) <- colnames(mat)
    return(as.vector(psd(plucked, phylo)[ ,1]))
  }
  calcSE <- function (x) {
    return (sd(x, na.rm = TRUE)/sqrt(length(x)))
  }
  # now for the actual calculations ...
  mat <- as.matrix(comm.data)
  if (metric == 'PD') {
    mat.phylo <- sapply(1:nrands, PDrarefy)
  } else if (metric == 'PSV') {
    mat.phylo <- sapply(1:nrands, PSVrarefy)
  } else if (metric == 'PSE') {
    mat.phylo <- sapply(1:nrands, PSErarefy)
  } else if (metric == 'PSR') {
    mat.phylo <- sapply(1:nrands, PSRrarefy)
  } else if (metric == 'PSC') {
    mat.phylo <- sapply(1:nrands, PSCrarefy)
  } else if (metric == 'PSD') {
    mat.phylo <- sapply(1:nrands, PSDrarefy)
  } else {
    stop('Unknown metric. Use: PD, PSV, PSE, PSR, PSC or PSD')
  }
  out.phylo <- apply(mat.phylo, 1, mean, na.rm = TRUE)
  out.se <- apply(mat.phylo, 1, calcSE)
  # certain metrics can produce NA
  out <- cbind(out.phylo, out.se)
  colnames(out) <- c(metric, paste0(metric, '.se'))
  return(out)
}

#phylodist (Author: J. Lim)
# Calculates phylogenetic nearest native neighbour distance given a list of alien species and native species
# Also calculates the mean phylogenetic distance between each alien species and all the native species on the list

phylodist <- function(natives, aliens, phy){
  #Diagnostic checks
  taxon <- c(natives, aliens)
  if(sum(taxon %in% phy$tip.label) != length(taxon)){
    print("Phylogeny does not contain all the taxa. PNND for missing taxa cannot be calculated and will be ignored.")
  }
  tip <- phy$tip.label[!phy$tip.label %in% taxon] 
  trimmedphy <- drop.tip(phy, tip = tip) #Removing taxa that are not found; phylogeny ONLY contains natives and aliens
  alien.incl <- aliens %in% trimmedphy$tip.label
  natives.incl <- natives %in% trimmedphy$tip.label
  excludedaliens <- aliens[! alien.incl]
  excludednatives <- natives[! natives.incl] # Output unmatched species for sanity checking
  trimmedaliens <- aliens[alien.incl] #New aliens list which now no longer contains aliens not represented in the tree
  trimmednatives <- natives[natives.incl]
  trimmedtaxon <- c(trimmednatives, trimmedaliens)
  #Results vectors
  pnnd <- NULL
  mpd <- NULL
  closest.nat <- NULL
  #Distance matrices
  phydist <- cophenetic(trimmedphy) #Generate phylogenetic distance matrix
  order <- match(trimmedtaxon, colnames(phydist)) #Find index of rows for taxa
  phydist <- phydist[order, order] #Ensuring that the rows are sorted in the same order as taxa
  exclude <- which(colnames(phydist) %in% trimmedaliens == TRUE) #Find rows for aliens
  temp <- phydist[-exclude, ,drop = FALSE] #Remove all non-natives rows from matrix
  for(i in trimmedaliens){
    alien.ind <- colnames(temp) %in% i
    alien.col <- temp[, alien.ind, drop = FALSE]
    closest.native.ind <- which.min(alien.col) #Find row containing minimum phylogenetic distance along column; shouldn't be any other non-natives, or itself among the rows
    closest.nat <- c(closest.nat, rownames(temp)[closest.native.ind]) #Extract closest native taxon name
    pnnd <- c(pnnd, temp[closest.native.ind, alien.ind]) #Extract phylogenetic distance
    mpd <- c(mpd, sum(alien.col) / dim(alien.col)[1])
  }
  results <- data.frame(pnnd = pnnd, mpd = mpd, taxon = trimmedaliens, closest.nat = closest.nat)
  output <- list(results = results, excluded.aliens = excludedaliens, excluded.natives = excludednatives)
}

#Creates a list of trait matrices
traitdistlist <- function(trait){
  traitlist <- names(trait)
  ntrait <- length(traitlist)
  traitdist <- vector("list", ntrait)
  names(traitdist) <- colnames(trait)
  for(j in 1:ntrait){
    trait.col <- names(trait) %in% traitlist[j]
    if(is.factor(trait[, trait.col]) == TRUE){
      levels(trait[, trait.col]) <- 1:length(levels(trait[, trait.col])) #Convert levels of trait into numerical values
      temp <- dist(trait[,j], method = "manhattan", upper = TRUE, diag = TRUE)
      temp <- as.matrix(temp)
      temp <- ifelse(temp > 0, 1, 0)
      rownames(temp) <- colnames(temp) <- rownames(trait)
      traitdist[[j]] <- temp
    } else {
      temp <- dist(trait[, trait.col], method = "manhattan", upper = TRUE, diag = TRUE)
      temp <- as.matrix(temp)
      rownames(temp) <- colnames(temp) <- rownames(trait)
      traitdist[[j]] <- temp
    }
  }
  return(traitdist)
}

# Calculates trait distance based on closest native
traitdist <- function(closest.native, aliens, traitdistlist){
  stopifnot(length(aliens) == length(closest.native)) #Error catching
  traitdiff <- NULL
  for(i in 1:length(aliens)){
    for(j in 1:length(traitdistlist)){
      alien.col <- which(rownames(traitdistlist[[j]]) %in% aliens[i] == TRUE)
      native.row <- which(rownames(traitdistlist[[j]]) %in% closest.native[i] == TRUE)
      traitdiff <- c(traitdiff, traitdistlist[[j]][native.row, alien.col])
    }
  }
  traitresults <- matrix(data = traitdiff, nrow = length(aliens), ncol = length(traitdistlist), byrow = TRUE)
  colnames(traitresults) <- paste(names(traitdistlist), ".diff", sep = "")
  traitresults <- as.data.frame(traitresults)
  traitresults$taxon <- aliens
  return(traitresults)
}

taxa2phylomatic <- function(taxa, output.dir, output.name, binom, storedtree){
  storedtree= "R20120829"
  informat = outformat = "newick"
  method = "phylomatic"
  url = 'http://phylodiversity.net/phylomatic/pmws'
  taxaformat = "slashpath"
  constraint.dir = file.path(output.dir, paste(output.name, "_tree.nwk", sep =""))
  unmatched.dir = file.path(output.dir, paste(output.name, "_unmatched.txt", sep = ""))
  # Phylomatic does not accept line breaks
  taxa <- str_replace_all(taxa, " ", "_")
  binom <- str_replace_all(binom, " ", "_")
  taxa <- paste(taxa, collapse = "\n")
  taxastring <- str_replace_all(str_replace_all(taxa, "/", "%2F"), "\n", "%0D")
  # Writing phylomatic query string
  phylomat.dat1 = paste("storedtree=", storedtree, "&", "informat=", informat, "&", "method=", method, "&", sep = "") 
  phylomat.dat2 = paste("taxaformat=", taxaformat, "&", "clean=true", "&", "outformat=", outformat, "&", sep = "")
  phylomat.dat3 = paste("taxa=", taxastring, sep = "")
  fullstr <- paste(phylomat.dat1, phylomat.dat2, phylomat.dat3, sep = "")
  # Parsing string to phylomatic using Curl
  print("Generating constraint tree using phylomatic...")
  temp <- postForm(uri = url, curl = getCurlHandle(), .opts = list(postfields = fullstr))
  # Cleaning up names (because tip labels may change in capitalization)
  for(i in 1:length(binom)){
    print(binom[i])
    temp <- gsub(pattern = binom[i], x = temp, replacement = binom[i], ignore.case = TRUE)
  }
  # Removing unmatched taxa
  matched <- gsub("\\[.+\\]", "", temp)
  # Export newick file
  print(paste("Writing constraint tree to ", constraint.dir))
  constraintnwk <-file(constraint.dir)
  writeLines(matched, constraintnwk)
  close(constraintnwk)
  # Exporting unmatched taxa to separte document
  m <- regexpr("\\[.+\\]", temp)
  unmatched <- regmatches(temp, m)
  print(unmatched)
  print(paste("Writing unmatched taxa to ", unmatched.dir))
  unmatchedfile <- file(unmatched.dir)
  writeLines(unmatched, unmatchedfile)
  close(unmatchedfile)
  #Defunct code. Does not handle large Request-URIs
  #urlquery <- paste(url, "?", phylomat.dat1, phylomat.dat2, phylomat.dat3, sep = "")
  #tt <- getURLContent(urlquery, curl = getCurlHandle())
  #phylomatic_nwk <- tt[[1]]
  #return(phylomatic_nwk)
  
  ## TO DO##
  #Argument for phylomatic tree to use
  #Argument to supply user trees
}

#newick and taxo2newick (Author: J. Lim)
#newick returns a list of taxa names into a polytomic newick tree
#taxo2newick returns a dataframe of taxa names into a newick tree,
#assuming monophyly in the penultimate taxonomic hierarchy (e.g., if data
#is on genera, then confamilial genera will be a monophyletic polytomy)
newick <- function(x){
  base <- x[1]
  if(length(x) > 1){
    for(i in 2:length(x)){
      base <- paste(base, x[i], sep = ", ")
    }
    return(paste("(", base, ")", sep = ""))
  } else {
    return(base)
  }
}

taxo2newick <-function(x, backbone = TRUE, backbonenewick, hightaxo, lowtaxo){
  taxo <- c("order", "class", "family")
  if(sum(match(tolower(names(x)), taxo, nomatch=0)) == 0){
    stop("Dataframe does not have taxonomic classifications")
  }else{
    hightaxoindex <- match(hightaxo, names(x))   #Find index of desired taxa
    hightaxolist <- as.vector(unique(x[, hightaxoindex]))   #Figure out what are the levels in the particular
    if(backbone == TRUE){
      base <- backbonenewick
    } else {
      base <- newick(hightaxolist)
    }
    for(i in 1:length(hightaxolist)){
      temp <- subset(x, x[,hightaxoindex] == hightaxolist[i]) #Subsets dataframe based on higher taxonomy to find the lower taxonomic levels
      lowtaxolist <- as.vector(unique(temp[,match(lowtaxo, names(temp))]))
      temp2 <- newick(lowtaxolist)
      base <- gsub(hightaxolist[i], temp2, base)
      print(base)
    }
  }
  return(base)
}

#osgf_to_en (Author: J. Lim)
#Finds centroids of British Ordinance Survey grid references and converts them to Eastings and Northings
#Accepts a vector of grid references and returns a data frame of eastings and northings 
#Works with both 2 (100x100 km), 4 (10x10 km), 6 (1x1 km) charactodigit grid reference numbers

osgr_to_en <- function(x, centroid = TRUE){

hectads <- c("HP", "HT", "HU", "HW", "HX", "HY", "HZ", "NA", "NB", "NC", "ND", "NF", "NG", "NH", "NJ", "NK", "NL", "NM", "NN", "NO", "NR", "NS", "NT", "NU", "NW", "NX", "NY", "NZ", "OV", "SC", "SD", "SE", "SH", "SJ", "SK", "SM", "SN", "SO", "SP", "SR", "SS", "ST", "SU", "SV", "SW", "SX", "SY", "SZ", "TA", "TF", "TG", "TL", "TM", "TQ", "TR", "TV")
X <- c(4, 3, 4, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 5, 6, 5, 6, 5, 6, 5)
Y <- c(12,11,11,10,10,10,10,9, 9, 9, 9, 8, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 4, 3, 3, 2, 2, 1, 1, 0)

res <- nchar(x[1])
eastings <- rep(NA, length(x))
northings <- rep(NA, length(x))
gridletter <- rep(NA, length(x))
x2 <- x1 <- rep(NA, length(x))
y2 <- y1 <- rep(NA, length(x))

if(res == 2){
	for(i in 1:length(x)){
		gridletter[i] <- substr(x[i], start = 1, stop = 2) #100 x 100 km scale
	}
}
if(res == 4){
	for(i in 1:length(x)){
		gridletter[i] <- substr(x[i], start = 1, stop = 2) #100 x 100 km scale
		x1[i] <- as.numeric(substr(x[i], start = 3, stop = 3)) #10 x 10 scale
		y1[i] <- as.numeric(substr(x[i], start = 4, stop = 4))
	}
} 
if(res == 6){
	for(i in 1:length(x)){
		gridletter[i] <- substr(x[i], start = 1, stop = 2) #100 x 100 km scale
		x1[i] <- as.numeric(substr(x[i], start = 3, stop = 3)) #10 x 10 scale
		x2[i] <- as.numeric(substr(x[i], start = 4, stop = 4))
		y1[i] <- as.numeric(substr(x[i], start = 5, stop = 5))
		y2[i] <- as.numeric(substr(x[i], start = 6, stop = 6))
	}
}

	for(i in 1:length(x)){
		if(res == 2){
			index <- match(gridletter[i], hectads)
			eastings[i] <- (X[index] * 100000)
			northings[i] <- (Y[index] * 100000)
		}
		if(res == 4){
			index <- match(gridletter[i], hectads)
			eastings[i] <- (X[index] * 100000) + (x1[i] * 10000)
			northings[i] <- (Y[index] * 100000) + (y1[i] * 10000)	
			} else {
			index <- match(gridletter[i], hectads)
			eastings[i] <- (X[index] * 100000) + (x1[i] * 10000) + (x2[i] * 1000) 
			northings[i] <- (Y[index] * 100000) + (y1[i] * 10000) + (x2[i] * 1000)
			}
	}
	if(centroid == TRUE){
		if(res == 2){
			eastings <- eastings + 50000
			northings <- northings + 50000
		}
		if(res == 4){
			eastings <- eastings + 5000
			northings <- northings + 5000
		}
		if(res == 6){
			eastings <- eastings + 500
			northings <- northings + 500
		}
	}
	data.frame(hectad = x, eastings, northings)
}

## RANDOMIZATION TEST FOR COLLESS Ic METRIC OF PHYLOGENETIC IMBALANCE
randI <- function(splist, phy, nrand){
  nsp <- length(splist)
  tip <- phy$tip.label[! phy$tip.label %in% splist]
  obsphy <- drop.tip(phy, tip = tip)
  obsphy <- as.treeshape(obsphy)
  obsI <- colless(obsphy, norm = "pda")
  randI <- NULL
  for(i in 1:nrand){
    randsp <- sample(phy$tip.label, size=nsp, replace=FALSE)
    tip <- phy$tip.label[! phy$tip.label %in% randsp]
    tempphy <- drop.tip(phy, tip)
    tempphy <- as.treeshape(tempphy)
    temp <- colless(tempphy, norm = "pda")
    randI <- c(randI, temp)
  }
  rank <- rank(c(obsI, randI), ties.method="max")[1]
  results <- data.frame(splocal = nsp, spreg = length(phy$tip.label), nrand = nrand, localcolless = obsI, p = rank/nrand)
  output <- list(results = results, randI = randI)
  return(output)
}
