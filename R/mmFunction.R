##### Function to conduct trial analysis ---


trial.analysis = function(df = df, Value = Value , Occ = Occ, Gen_no = Gen_no, 
                          Rep = Rep, Sub_block = Sub_block) {
  if (nlevels(df$Occ)>1) {
    lmm.1 = lmer(Value ~ Occ + 
                   (1|Gen_no) + (1|Gen_no:Occ) + 
                   (1|Rep) + (1|Occ:Rep:Sub_block),
                 data = df)
  } else {
    lmm.2 = lmer(Value ~ 1 + 
                   (1|Gen_no)+ (1|Rep) + (1|Rep:Sub_block),
                 data = df)
  }
}
