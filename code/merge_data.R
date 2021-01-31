# sessions 1-19
df1 = read.csv("../data/1_19.csv",header=T)
# sessions 20-38
df2 = read.csv("../data/20_38.csv",header=T)
# sessions 39-59
df3 = read.csv("../data/39_59.csv",header=T)
# sessions 60-77
df4 = read.csv("../data/60_77.csv",header=T)

# modularities for all sessions
mod_df = read.csv("../data/modularities.csv",header=T)

# full network dataframe
full_net_df = cbind(rbind(df1,df2,df3,df4),mod_df)

# behavioral data
b_df = read.csv("../data/behavioral_scores.csv",header=T)
b_df = b_df[-c(1,21,52,58),-1]

# joined dataframe
full_df = cbind(full_net_df,b_df)
full_df = full_df[,-1]

#Transform the label cluster coeff and Q with Fisher_to_Z 
r=full_df[,2:3]
full_df[,2:3] = 0.5 * log((r)/(1-r))

