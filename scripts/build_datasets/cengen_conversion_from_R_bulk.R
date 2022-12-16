cengen_sc_1_bulk <-
  local({
    load(paste(datafolder,"cengen_sc_1.rda",sep=""))

    new_columns <- cbind(
      AWC = apply(cengen_sc_1[,c("AWC_ON", "AWC_OFF")], 1, mean),
      DA  = apply(cengen_sc_1[,c("DA", "DA9")], 1, mean),
      DB  = apply(cengen_sc_1[,c("DB", "DB01")], 1, mean),
      IL2 = apply(cengen_sc_1[,c("IL2_DV", "IL2_LR")], 1, mean),
      RMD = apply(cengen_sc_1[,c("RMD_DV", "RMD_LR")], 1, mean),
      RME = apply(cengen_sc_1[,c("RME_DV", "RME_LR")], 1, mean),
      VA  = apply(cengen_sc_1[,c("VA", "VA12")], 1, mean),
      VB  = apply(cengen_sc_1[,c("VB", "VB01", "VB02")], 1, mean),
      VC  = apply(cengen_sc_1[,c("VC", "VC_4_5")], 1, mean),
      VD  = cengen_sc_1[,"VD_DD"],
      DD  = cengen_sc_1[,"VD_DD"])

    redundant_columns <- c("AWC_ON","AWC_OFF","DA","DA9","DB","DB01","IL2_DV","IL2_LR",
                           "RMD_DV","RMD_LR","RME_DV","RME_LR","VA","VA12","VB","VB01","VB02",
                           "VC","VC_4_5","VD_DD")

    res <- cbind(
      cengen_sc_1[ , - which(colnames(cengen_sc_1) %in% redundant_columns)],
      new_columns
    )

    # alphabetic order
    res[,order(colnames(res))]
  })
cengen_sc_2_bulk <-
  local({
    load(paste(datafolder,"cengen_sc_2.rda",sep=""))

    new_columns <- cbind(
      AWC = apply(cengen_sc_2[,c("AWC_ON", "AWC_OFF")], 1, mean),
      DA  = apply(cengen_sc_2[,c("DA", "DA9")], 1, mean),
      DB  = apply(cengen_sc_2[,c("DB", "DB01")], 1, mean),
      IL2 = apply(cengen_sc_2[,c("IL2_DV", "IL2_LR")], 1, mean),
      RMD = apply(cengen_sc_2[,c("RMD_DV", "RMD_LR")], 1, mean),
      RME = apply(cengen_sc_2[,c("RME_DV", "RME_LR")], 1, mean),
      VA  = apply(cengen_sc_2[,c("VA", "VA12")], 1, mean),
      VB  = apply(cengen_sc_2[,c("VB", "VB01", "VB02")], 1, mean),
      VC  = apply(cengen_sc_2[,c("VC", "VC_4_5")], 1, mean),
      VD  = cengen_sc_2[,"VD_DD"],
      DD  = cengen_sc_2[,"VD_DD"])

    redundant_columns <- c("AWC_ON","AWC_OFF","DA","DA9","DB","DB01","IL2_DV","IL2_LR",
                           "RMD_DV","RMD_LR","RME_DV","RME_LR","VA","VA12","VB","VB01","VB02",
                           "VC","VC_4_5","VD_DD")

    res <- cbind(
      cengen_sc_2[ , - which(colnames(cengen_sc_2) %in% redundant_columns)],
      new_columns
    )

    # alphabetic order
    res[,order(colnames(res))]
  })
cengen_sc_3_bulk <-
  local({
    load(paste(datafolder,"cengen_sc_3.rda",sep=""))

    new_columns <- cbind(
      AWC = apply(cengen_sc_3[,c("AWC_ON", "AWC_OFF")], 1, mean),
      DA  = apply(cengen_sc_3[,c("DA", "DA9")], 1, mean),
      DB  = apply(cengen_sc_3[,c("DB", "DB01")], 1, mean),
      IL2 = apply(cengen_sc_3[,c("IL2_DV", "IL2_LR")], 1, mean),
      RMD = apply(cengen_sc_3[,c("RMD_DV", "RMD_LR")], 1, mean),
      RME = apply(cengen_sc_3[,c("RME_DV", "RME_LR")], 1, mean),
      VA  = apply(cengen_sc_3[,c("VA", "VA12")], 1, mean),
      VB  = apply(cengen_sc_3[,c("VB", "VB01", "VB02")], 1, mean),
      VC  = apply(cengen_sc_3[,c("VC", "VC_4_5")], 1, mean),
      VD  = cengen_sc_3[,"VD_DD"],
      DD  = cengen_sc_3[,"VD_DD"])

    redundant_columns <- c("AWC_ON","AWC_OFF","DA","DA9","DB","DB01","IL2_DV","IL2_LR",
                           "RMD_DV","RMD_LR","RME_DV","RME_LR","VA","VA12","VB","VB01","VB02",
                           "VC","VC_4_5","VD_DD")

    res <- cbind(
      cengen_sc_3[ , - which(colnames(cengen_sc_3) %in% redundant_columns)],
      new_columns
    )

    # alphabetic order
    res[,order(colnames(res))]
  })
  
cengen_sc_4_bulk <-
  local({
    load(paste(datafolder,"cengen_sc_4.rda",sep=""))

    new_columns <- cbind(
      AWC = apply(cengen_sc_4[,c("AWC_ON", "AWC_OFF")], 1, mean),
      DA  = apply(cengen_sc_4[,c("DA", "DA9")], 1, mean),
      DB  = apply(cengen_sc_4[,c("DB", "DB01")], 1, mean),
      IL2 = apply(cengen_sc_4[,c("IL2_DV", "IL2_LR")], 1, mean),
      RMD = apply(cengen_sc_4[,c("RMD_DV", "RMD_LR")], 1, mean),
      RME = apply(cengen_sc_4[,c("RME_DV", "RME_LR")], 1, mean),
      VA  = apply(cengen_sc_4[,c("VA", "VA12")], 1, mean),
      VB  = apply(cengen_sc_4[,c("VB", "VB01", "VB02")], 1, mean),
      VC  = apply(cengen_sc_4[,c("VC", "VC_4_5")], 1, mean),
      VD  = cengen_sc_4[,"VD_DD"],
      DD  = cengen_sc_4[,"VD_DD"])

    redundant_columns <- c("AWC_ON","AWC_OFF","DA","DA9","DB","DB01","IL2_DV","IL2_LR",
                           "RMD_DV","RMD_LR","RME_DV","RME_LR","VA","VA12","VB","VB01","VB02",
                           "VC","VC_4_5","VD_DD")

    res <- cbind(
      cengen_sc_4[ , - which(colnames(cengen_sc_4) %in% redundant_columns)],
      new_columns
    )

    # alphabetic order
    res[,order(colnames(res))]
  })
