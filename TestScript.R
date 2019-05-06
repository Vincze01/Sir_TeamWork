

newBorn = list()
timeToVaccination = 5

for(i in 1:21)
{
  oldestItem = NULL
  if(length(newBorn) >= timeToVaccination)
  {
    oldestItem = unlist(newBorn[1])
    newBorn[1] = NULL
  }
  newBorn = c(newBorn, i)
  
  for (index in 1:length(newBorn))
  {
    newBorn[index] = unlist(newBorn[index]) * 2
  }
  
  if(is.null(oldestItem) == FALSE) #Csak ha mar van ilyen elem akkor csinalsz vele barmit
  {
    #itt
  }
  
  for (item in newBorn) {
    print(item)
  }
  print("###########################")  
}
