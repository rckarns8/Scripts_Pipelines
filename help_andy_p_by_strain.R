library(dplyr)
library(tidyr)
library(hash)

test_var <- c('nh4', 'nh4c', 'mox')
p_values <- hash()
max_length <- 0

grouped_data_for_table <- aoa %>%
  distinct(treat, time) %>%
  arrange(treat, time)


for (i in test_var) {
  p_values[[i]] = list()
  #Need to convert data in to strain, treat, time with nh4 data
  #combined in to a list.
  #Note: In the final grouped_data variable, the column "value" represents
  #nh4 data.
  grouped_data <- aoa %>% 
    #To change data that is used in t test, swap the nh4 variable below
    #for your new column name
    select(treat, strain, time, i) %>%
    gather(key = variable, value = value, -treat, -strain, -time)  %>%
    group_by(treat, strain, time) %>% 
    summarise(value = list(value))
  
  #Since we now have grouped data, separate
  #time 1 from time 2. At this point strain_d[j] = strain_n[j] for treat and strain
  strain_d <- subset(grouped_data, strain == 'D')
  strain_d <- subset(strain_d, treat != 'HH')
  strain_d <- subset(strain_d, treat != 'L')
  strain_n <- subset(grouped_data, strain == 'N')
  
  #Set length of loop to avoid duplicates in final table
  max_length <- nrow(strain_d)
  
  #For loop
  for (j in seq(nrow(strain_d))){
    #Grab i (variable from test_var) values from strain_d and strain_n Remove NA values and then cast list as numeric
    temp_1 <- as.numeric(lapply(strain_d$value[[j]], function(x) x[!is.na(x)]))
    temp_2 <- as.numeric(lapply(strain_n$value[[j]], function(x) x[!is.na(x)]))
    p_values[[i]] = c(p_values[[i]], list(t.test(temp_1, temp_2, paired = TRUE)$p.value))
  }
}



#Create empty data frame for treat, strain, and p_value
final_p_values_strain <- data.frame(treat=character(), strain=character(), nh4_p_value=double(), nh4c_p_value=double(), mox_p_value=double())

#Set our dataframe (final_p_values) to a new data frame that is
#the current final_p_values data frame merged with a "new" dataframe which 
#represents our new row. This row will contain:
#strain: time_1[j] (which is also time_2[j])
#treat: time_1[j] (which is also time_2[j])
#p_value: result from the t test using time_1[j] and time_2[j] nh4 values
for (j in seq(max_length)){
  final_p_values_strain <- rbind(final_p_values_strain, data.frame(
    treat = grouped_data_for_table$treat[[j+1]],
    time = grouped_data_for_table$time[[j+1]],
    nh4_p_value=as.numeric(p_values[["nh4"]][j]),
    nh4c_p_value=as.numeric(p_values[["nh4c"]][j]),
    mox_p_value=as.numeric(p_values[["mox"]][j])
  ))
}
