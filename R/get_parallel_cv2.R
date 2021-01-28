get_parallel_cv2 <- function(X,
                                Y,
                                lambdas,
                                non_zeros,
                                label,
                                penalization,
                                max_iterations,
                                tolerance,
                                modes2){
  #non-parallel x-fold Cross validation for each alpha, for example alpha = 0, 0.1, ... , 0.9, 1.0
  #
  #input:     X             n*p matrix          - independent variables
  #           Y             n*q matrix          - dependent variables
  #           lambda        integer             - factor for Ridge penalty
  #           labels         vector with n_subset unique elements with length of n - for crossvalidation
  #
  #output:    abs_cor       vector of integers  - returns the summed abs correlation
  #           stime         int/time            - return algorithms running time

  #magic variables for testing
  # lambdas <- lambda
  # non_zeros <- nonzero
  # X = X.sampled
  # Y = Y.sampled
  # tolerance = tol
  # max_iterations = maxiter
  #********************

  nr_subsets      <-    length(unique(label))
  nmodes          <-    modes2
  abs_cors        <-    c()
  iterations_m    <-    c()

  length(lambdas)
  length(non_zeros)

  #on windows
  # cl = makeCluster(no_cores)
  # clusterExport(cl, c("X","label",
  #                     "Y", "splspm",
  #                     "check_args", "check_data",
  #                     "is_not_tabular", "lacks_rownames", "lacks_colnames",
  #                     "check_path", "is_not_matrix", "is_square_matrix"),
  #               envir=environment())

  sub_abs_cor       <- c()
  sub_results       <- c()
  sub_iterations_m  <- c()

  closeAllConnections()
  no_cores <- detectCores() - 1
  cl = makeCluster(no_cores, type = "FORK")

  ## RDA objective function
  if(nmodes[1]!="B" || nmodes[2]!="B"){
    #measure time
    stime <- system.time({
    par_CVs_correlations <- parSapply(cl, 1:length(non_zeros),
                                  function(nz){
                                    par_sub_abs_cor <- sapply(1:length(lambdas), function(l){

                                      nz_correlations <- sapply(1:nr_subsets,
                                                                    function(i){
                                                                      seed_nr <- set.seed(i, kind = "L'Ecuyer-CMRG")
                                                                      X.train   <- X[label!=i,]
                                                                      X.test    <- X[label==i,]
                                                                      X.train   =   X.train[,apply(X.test,2,sd)>0]
                                                                      X.test    =   X.test[,apply(X.test,2,sd)>0]
                                                                      X.test    =   X.test[,apply(X.train,2,sd)>0]
                                                                      X.train   =   X.train[,apply(X.train,2,sd)>0]
                                                                      Y.train   <- Y[label!=i,]
                                                                      Y.test    <- Y[label==i,]
                                                                      Y.train   =   Y.train[,apply(Y.test,2,sd)>0]
                                                                      Y.test    =   Y.test[,apply(Y.test,2,sd)>0]
                                                                      Y.test    =   Y.test[,apply(Y.train,2,sd)>0]
                                                                      Y.train   =   Y.train[,apply(Y.train,2,sd)>0]
                                                                      ndata_sets <- cbind(X.train,Y.train)
                                                                      nEXPL_X = c(0,0)
                                                                      nRESP_Y = c(1,0)
                                                                      npath_matrix = rbind(nEXPL_X, nRESP_Y)
                                                                      nblocks = list(1:dim(X.train)[2], dim(X.train)[2]+1:dim(Y.train)[2])
                                                                      sub_results <- splspm(ndata_sets,npath_matrix, nblocks,
                                                                                            nmodes,scheme="path",penalization = penalization,
                                                                                            nonzero = non_zeros[nz],lambda = lambdas[l],
                                                                                            maxiter = max_iterations,tol = tolerance,
                                                                                            warning_non_convergence = FALSE)
                                                                      ALPHA <- t(X.train)%*% sub_results$scores[,1]%*%solve(t( sub_results$scores[,1])%*%sub_results$scores[,1])
                                                                      XI.test = scale(X.test) %*% ALPHA
                                                                      sub_abs_cor <- sum((abs(cor(XI.test,Y.test))))/dim(Y.train)[2]
                                                                      c(sub_abs_cor = sub_abs_cor)
                                                                      })
                                      c(nz_correlations = nz_correlations)
                                      })#end of non_zeros sapply
                                  c(par_sub_abs_cor = par_sub_abs_cor)
                                  }) #end of lambda sapply
    })[3]#end of measure time

    results.par <- do.call(rbind,as.list(par_CVs_correlations))
    results.par <- lapply(c(1:(length(non_zeros)*length(lambdas))),function(x) results.par[(((nr_subsets*x)-nr_subsets)+1):(nr_subsets*x),])
    results.par <- do.call(cbind,results.par)
    results.par

  }#END OF RDA style
  else{
    #START OF CCA STYLE

    #measure time
    stime <- system.time({
      par_CVs_correlations <- parSapply(cl, 1:length(non_zeros),
                                        function(nz){
                                          par_sub_abs_cor <- sapply(1:length(lambdas), function(l){
                                            nz_correlations <- sapply(1:nr_subsets,
                                                                      function(i){
                                                                        seed_nr <- set.seed(i, kind = "L'Ecuyer-CMRG")
                                                                        X.train   <- X[label!=i,]
                                                                        X.test    <- X[label==i,]
                                                                        X.train   =   X.train[,apply(X.test,2,sd)>0]
                                                                        X.test    =   X.test[,apply(X.test,2,sd)>0]
                                                                        X.test    =   X.test[,apply(X.train,2,sd)>0]
                                                                        X.train   =   X.train[,apply(X.train,2,sd)>0]
                                                                        Y.train   <- Y[label!=i,]
                                                                        Y.test    <- Y[label==i,]
                                                                        Y.train   =   Y.train[,apply(Y.test,2,sd)>0]
                                                                        Y.test    =   Y.test[,apply(Y.test,2,sd)>0]
                                                                        Y.test    =   Y.test[,apply(Y.train,2,sd)>0]
                                                                        Y.train   =   Y.train[,apply(Y.train,2,sd)>0]
                                                                        ndata_sets <- cbind(X.train,Y.train)
                                                                        nEXPL_X = c(0,0)
                                                                        nRESP_Y = c(1,0)
                                                                        npath_matrix = rbind(nEXPL_X, nRESP_Y)
                                                                        nblocks = list(1:dim(X.train)[2], dim(X.train)[2]+1:dim(Y.train)[2])
                                                                        sub_results <- splspm(ndata_sets,npath_matrix,nblocks,
                                                                                              nmodes,scheme="path",penalization = penalization,
                                                                                              nonzero = non_zeros[nz],lambda = lambdas[l],
                                                                                              maxiter = max_iterations,tol = tolerance,
                                                                                              warning_non_convergence = FALSE,cross_validate = FALSE)
                                                                        ALPHA <- t(X.train)%*% sub_results$scores[,1]%*%
                                                                          solve(t( sub_results$scores[,1])%*% sub_results$scores[,1])
                                                                        XI.test = scale(X.test) %*% ALPHA
                                                                        BETA <- t(Y.train)%*% sub_results$scores[,2]%*%
                                                                          solve(t( sub_results$scores[,2])%*%sub_results$scores[,2])
                                                                        ETA.test = scale(Y.test) %*% BETA
                                                                        sub_abs_cor <- sum((abs(cor(XI.test,ETA.test))))
                                                                        c(sub_abs_cor = sub_abs_cor)})
                                            c(nz_correlations = nz_correlations)})#end of non_zeros sapply
                                          c(par_sub_abs_cor = par_sub_abs_cor)
                                        }) #end of lambda sapply
    })[3]#end of measure time

    results.par <- do.call(rbind,as.list(par_CVs_correlations))
    results.par <- lapply(c(1:(length(non_zeros)*length(lambdas))),function(x) results.par[(((nr_subsets*x)-nr_subsets)+1):(nr_subsets*x),])
    results.par <- do.call(cbind,results.par)
    results.par

    }

  abs_cors <- results.par
  #Figure out lambdas and non-zeros columns in results
  labels_non_zeros  <- rep(non_zeros, dim(abs_cors)[2]/length(non_zeros))
  labels_non_zeros

  labels_lambdas    <- rep(lambdas, each=dim(abs_cors)[2]/length(lambdas))
  labels_lambdas
  length(labels_lambdas)

  all_abs_cors  <- rbind(labels_lambdas, labels_non_zeros, abs_cors)


  mean_abs_cors <- c()

  for (i in 1:length(labels_lambdas)){

    sub_result    <-  (c(labels_lambdas[i], labels_non_zeros[i], mean(abs_cors[,i], na.rm = TRUE)))
    mean_abs_cors <- rbind(mean_abs_cors, sub_result)

  }
  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Ridge Penalty",
                                 "Number of Nonzeros", "mean_abs_cors")
  mean_abs_cors


  #plot(mean_abs_cors[,1],mean_abs_cors[,3], pch=19, col=mean_abs_cors[,2])
  #text(mean_abs_cors[,1],mean_abs_cors[,3], labels=mean_abs_cors[,2], cex= 1.1, pos=2, pch=19, col=mean_abs_cors[,2])


  # #*********************************
  #
  # plot2 <-
  #   ggplot(data=data.frame(mean_abs_cors),
  #          aes(x = factor(Lambda), y = mean_abs_cors,
  #              group = factor(Number.of.Nonzeros),
  #              shape = factor(Number.of.Nonzeros),
  #              color = factor(Number.of.Nonzeros)))+
  #   geom_line() +
  #   geom_point() +
  #   scale_x_discrete("Lambda") +
  #   scale_y_continuous("Mean absolute correlation") +
  #   facet_grid(.~Number.of.Nonzeros )
  #
  # #*********************************

  print("Elapsed time")
  print(stime)
  closeAllConnections()
  #Return section**********************
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime,
                           nmodes = nmodes
                           # plot2 = plot2,
                           #iterations_m = iterations_m
  )

  result
}
