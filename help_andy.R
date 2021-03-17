library(dplyr)
library(tidyr)
library(hash)

#Need to convert data in to strain, treat, time with nh4 data
#combined in to a list.
#Note: In the final grouped_data variable, the column "value" represents
#nh4 data.
grouped_data <- aoa %>% 
  #To change data that is used in t test, swap the nh4 variable below
  #for your new column name
  select(treat, strain, time, test_var) %>%
  gather(key = variable, value = value, -treat, -strain, -time)  %>%
  group_by(treat, strain, time) %>% 
  summarise(value = list(value))

#Since we now have grouped data, separate
#time 1 from time 2. At this point time_1[j] = time_2[j] for treat and strain
time_1 <- subset(grouped_data, time == 1)
time_2 <- subset(grouped_data, time == 2)

strain_d <- subset(grouped_data, strain == 'D')
strain_n <- subset(grouped_data, strain == 'N')

#Create empty data frame for treat, strain, and p_value
testresults <- data.frame(treat=character(), strain=character(), p_value=double())
#For loop
for (j in seq(nrow(time_1))){
  #Grab nh4 values from time_1 and time_2. Remove NA values and then cast list as numeric
  temp_1 <- as.numeric(lapply(time_1$value[[j]], function(x) x[!is.na(x)]))
  temp_2 <- as.numeric(lapply(time_2$value[[j]], function(x) x[!is.na(x)]))

  #Set our dataframe (testresults) to a new data frame that is
  #the current testresults data frame merged with a "new" dataframe which 
  #represents our new row. This row will contain:
  #strain: time_1[j] (which is also time_2[j])
  #treat: time_1[j] (which is also time_2[j])
  #p_value: result from the t test using time_1[j] and time_2[j] nh4 values
  testresults <- rbind(testresults, data.frame(
    treat = time_1$treat[[j]],
    strain = time_1$strain[[j]],
    p_value = t.test(temp_1, temp_2, paired = TRUE)$p.value
  ))
}
