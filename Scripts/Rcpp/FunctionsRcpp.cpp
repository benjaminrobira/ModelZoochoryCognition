#include <Rcpp.h>
#include <math.h> // Allows the use of mathematical functions
//Note: for fmod, because doing the modulo keeping sign (i.e. respecting mathematic def), I always add cycleLength, the maximum possible length, in order to avoid keeping the negative value instead of having the related positive one)
#include <cmath> // For abs and NA use
#include <iostream> // Allows to display input and output in the console
#include <fstream> // Allows to add/extract input/output to a file
#include <stdio.h>
#include <string> // Allow to use string
#include <string.h>// To use memcpy

using namespace Rcpp;

//----- General

struct add_multiple {
  int incr;
  int count;
  add_multiple(int incr)
    : incr(incr), count(0)
  {}
  inline int operator()(int d) {
    return d + incr * count++;
  }
};

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_seq(double from_, double to_, double by_ = 1.0) {
  
  if(from_ > to_){
    int transitory = from_;
    from_ = to_;
    to_ = transitory;
    by_ = abs(by_);
    
    int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
    int from = adjust * from_;
    int to = adjust * to_;
    int by = adjust * by_;
    
    std::size_t n = ((to - from) / by) + 1;
    Rcpp::IntegerVector res = Rcpp::rep(from, n);
    add_multiple ftor(by);
    
    std::transform(res.begin(), res.end(), res.begin(), ftor);
    return rev(Rcpp::NumericVector(res) / adjust);
  }else{
    int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
    int from = adjust * from_;
    int to = adjust * to_;
    int by = adjust * by_;
    
    std::size_t n = ((to - from) / by) + 1;
    Rcpp::IntegerVector res = Rcpp::rep(from, n);
    add_multiple ftor(by);
    
    std::transform(res.begin(), res.end(), res.begin(), ftor);
    return Rcpp::NumericVector(res) / adjust;
  }
}

//Find which max or min, available here: https://stackoverflow.com/questions/59055902/find-index-of-all-max-min-values-in-vector-in-rcpp

// [[Rcpp::export]]
IntegerVector which_maxRcpp(NumericVector v) {
  double current_max = v[0];
  int n = v.size();
  std::vector< int > res;
  res.push_back( 0 );
  int i;
  
  for( i = 1; i < n; ++i) {
    double x = v[i];
    if( x > current_max ) {
      res.clear();
      current_max = x;
      res.push_back( i );
    } else if ( x == current_max ) {
      res.push_back( i );
    }
  }
  Rcpp::IntegerVector iv( res.begin(), res.end() );
  //Return vector of positions of the maximum value in the vector
  return iv;
}

// [[Rcpp::export]]
IntegerVector which_minRcpp(NumericVector v) {
  double current_min = v[0];
  int n = v.size();
  std::vector< int > res;
  res.push_back( 0 );
  int i;
  
  for( i = 1; i < n; ++i) {
    double x = v[i];
    if( x < current_min ) {
      res.clear();
      current_min = x;
      res.push_back( i );
    } else if ( x == current_min ) {
      res.push_back( i );
    }
  }
  Rcpp::IntegerVector iv( res.begin(), res.end() );
  
  //Return vector of positions of the minimum value in the vector
  return iv;
}

// [[Rcpp::export]]
IntegerVector which_naRcpp(NumericVector v) {
  int n = v.size();
  std::vector< int > res;
  int i;
  
  for( i = 1; i < n; ++i) {
    double x = v[i];
    if( NumericVector::is_na(x) ) {
      res.push_back( i );
    }
  }
  Rcpp::IntegerVector iv( res.begin(), res.end() );
  
  //Return vector of positions of the NAs in the vector
  return iv;
}


//Modulo function for vector
// [[Rcpp::export]]
NumericVector moduloVector(
    NumericVector inputVector,
    double moduloValue
){
  NumericVector outputVector = inputVector;
  for(int count = 0; count < inputVector.size(); count++){
    outputVector[count] = fmod(outputVector[count] + 10*moduloValue, moduloValue);
  }
  
  //Return a vector of the modulos
  return outputVector;
}

//--------------------
// Calculate travelled distance
//--------------------

// This function calculates the beeline distance travelled by the agent

// [[Rcpp::export]]
double distanceTravelledCalculation(
    NumericVector previousCoordinates,
    NumericVector newCoordinates){
  double outputDistance = sqrt(
    pow(newCoordinates[0] - previousCoordinates[0], 2) +
      pow(newCoordinates[1] - previousCoordinates[1], 2)
  );
  
  //Return a value (double) indicating the distance travelled during the bout
  return(outputDistance);
}

//--------------------
// Get coordinates changing the referential so that travelled path = x axis
//--------------------

//This function rotates the coordinate axes so that the x axis is now the travelled beeline

// [[Rcpp::export]]
NumericMatrix changeReferentialCoordinatesAlongTravel(
    NumericVector departureLocation,
    NumericVector arrivalLocation,
    NumericMatrix treeLoc
){
  
  double distanceToTarget = distanceTravelledCalculation(
    departureLocation,
    arrivalLocation
  );
  
  NumericMatrix treeLocReferentialAgent(treeLoc(_,0).size(), 2);
  treeLocReferentialAgent(_,0) = treeLoc(_,0);
  treeLocReferentialAgent(_,1) = treeLoc(_,1);
  
  //I will apply roation to change referential so as to have classical orthogonal x and y axis to determine whether trees were seen by the agent
  //Change coordinates to have them in referential; Reminder: to do so, geometry formula can be found in wikipedia:
  //x coord
  treeLocReferentialAgent(_,0) =
    (treeLoc(_,0) - departureLocation[0]) *
    (arrivalLocation[0] - departureLocation[0])/distanceToTarget +
    (treeLoc(_,1) - departureLocation[1]) *
    (arrivalLocation[1] - departureLocation[1])/distanceToTarget;
  
  //y coord
  treeLocReferentialAgent(_,1) =
    -(treeLoc(_,0) - departureLocation[0]) *
    (arrivalLocation[1] - departureLocation[1])/distanceToTarget +
    (treeLoc(_,1) - departureLocation[1]) *
    (arrivalLocation[0] - departureLocation[0])/distanceToTarget;
  
  //Return two-column matrix about the coordinates in the "agent" referential
  return(treeLocReferentialAgent);
}   

//------------------------------------------------------------------------------------------------------------------------------------------------------------------


//----- Environment-related

//------------------
// Assign tree location
//------------------

//This function distributes the tree within the map. Distribution can be hommogeneous or heterogeneous.
//Heterogeneity is done by choosing some cluster center and distribution in a gaussian centered on this location (and the variance is parameterised with the cluster spread argument).

// [[Rcpp::export]]
NumericMatrix distributionTree(
    int numberTrees,
    double lowerBorder,//lower border of the map
    double upperBorder,//upper border of the map
    bool homogeneousDistribution = true,
    int treeClusterNumber = 0,
    int treeClusterSpread = 0
){
  
  NumericMatrix outputMatrix(numberTrees, 2);
  if(homogeneousDistribution){// Homogeneous distribution
    outputMatrix(_, 0) = runif(numberTrees,lowerBorder,upperBorder);
    outputMatrix(_, 1) = runif(numberTrees,lowerBorder,upperBorder);
  }else{// Heterogeneous distribution
    //Number of trees per cluster
    int numberTreePerCluster(numberTrees/treeClusterNumber);//Normally because both numerator and denominator are integers, it should also give an integer
    for(int count1 = 0; count1 < treeClusterNumber; count1++){
      //Calculate coordinates of each patch barycentre (homogeneous distribution of patches)
      double coordXPatch(R::runif(lowerBorder,upperBorder));//R space for having a scalar and not vector
      double coordYPatch(R::runif(lowerBorder,upperBorder));
      for(int count2 = 0; count2 < numberTreePerCluster; count2++){
        //Distribute trees around these patch barycentres
        outputMatrix(count1*numberTreePerCluster + count2, 0) = R::rnorm(coordXPatch, treeClusterSpread);
        //Ensure x is within map
        while(outputMatrix(count1*numberTreePerCluster + count2, 0) < lowerBorder || outputMatrix(count1*numberTreePerCluster + count2, 0) > upperBorder){
          outputMatrix(count1*numberTreePerCluster + count2, 0) = R::rnorm(coordXPatch, treeClusterSpread);
        }
        outputMatrix(count1*numberTreePerCluster + count2, 1) = R::rnorm(coordYPatch, treeClusterSpread);
        //Ensure y is within map
        while(outputMatrix(count1*numberTreePerCluster + count2, 1) < lowerBorder || outputMatrix(count1*numberTreePerCluster + count2, 1) > upperBorder){
          outputMatrix(count1*numberTreePerCluster + count2, 1) = R::rnorm(coordYPatch, treeClusterSpread);
        }
      }
      //Complete the number of trees that were not attributed to a cluster: homogeneous distribution for them
      int toCompleteTreeNumber(numberTrees - numberTreePerCluster*treeClusterNumber);
      for(int count3 = 0; count3 < toCompleteTreeNumber; count3++){
        outputMatrix((treeClusterNumber-1)*numberTreePerCluster + numberTreePerCluster + count3, 0) = R::runif(lowerBorder,upperBorder);
        outputMatrix((treeClusterNumber-1)*numberTreePerCluster + numberTreePerCluster + count3, 1) = R::runif(lowerBorder,upperBorder);
      }
    }
  }
  
  //Return a matrix with the x (first column) and y (second column) coordinates of tree position.
  return(outputMatrix);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------
// Assign the fruiting dates
//----------------------

//This function randomly (homogeneous distribution along time) assigns the fruiting dates accounting for circularity.

// [[Rcpp::export]]
NumericMatrix dateTree(
    int numberTrees,
    int cycleLength,
    int fruitingLength
){
  NumericMatrix outputMatrix(numberTrees, 2);
  outputMatrix(_,0) = runif(numberTrees, 0, cycleLength);
  outputMatrix(_,1) = outputMatrix(_,0) + fruitingLength;
  
  //Return matrix with start (first column) and end (second column) dates of fruiting
  return(outputMatrix);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Calculate food available at tree
//--------------------

//This function calculates the amount of food available based on what has been already eaten and the phenology.

// [[Rcpp::export]]
NumericVector foodTree(//Triangular distribution of food quantity
    double currentTime,
    NumericVector startTimeVector,//vector of fruiting start time for each tree
    NumericVector endTimeVector,//vector of fruiting end time for each tree
    double lengthOfCycle,
    NumericVector foodConsumedVector,
    NumericVector maximumFoodToYieldVector,//vector maximum food to be reached for a given tree (if saturation)
    bool linear = true//should the food linearily increase/decreae (as a triangle), with fruit rotting, or be constant?
){
  NumericVector outputVector(startTimeVector.size());
  NumericVector currentTimeVector =  rep(fmod(currentTime + 10*lengthOfCycle, lengthOfCycle), startTimeVector.size());
  
  //Translate so 0 matches start Time
  currentTimeVector = currentTimeVector - startTimeVector;
  endTimeVector = endTimeVector - startTimeVector;
  startTimeVector = startTimeVector - startTimeVector;
  
  // //Do the modulo (actually only necessary for timer since there has been a translation)
  currentTimeVector = moduloVector(currentTimeVector, lengthOfCycle);
  endTimeVector = moduloVector(endTimeVector, lengthOfCycle);
  startTimeVector = moduloVector(startTimeVector, lengthOfCycle);
  
  // Assign the food if no eating
  if(linear){
  outputVector = 
    ifelse(
      ((currentTimeVector > startTimeVector) & (currentTimeVector < endTimeVector)),
      ifelse(
        currentTimeVector <= startTimeVector + (endTimeVector - startTimeVector)/2,
                           maximumFoodToYieldVector * (currentTimeVector - startTimeVector)/(endTimeVector - startTimeVector) - foodConsumedVector,
                           maximumFoodToYieldVector - maximumFoodToYieldVector * (currentTimeVector - startTimeVector)/(endTimeVector - startTimeVector) - foodConsumedVector
      ),
      0.0 //No food                   
    );
  }else{
    outputVector = ifelse(
      ((currentTimeVector > startTimeVector) & (currentTimeVector < endTimeVector)),
      maximumFoodToYieldVector - foodConsumedVector,
      0.0 //No food                   
    );
  }
  
  //Readjust if negative values due to eating
  outputVector = 
    ifelse(
      outputVector < 0.0,
                     0.0,
                     outputVector
    );
  
  //Return a numeric vector with current food quantity yielded at each patch
  return(outputVector);
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// New coordinates for Dispersal
//--------------------

//This function assigns a new location along the way (buffer or not).

// [[Rcpp::export]]
NumericVector newCoordinatesAfterDispersal(
    NumericVector previousCoordinatesAgent,
    NumericVector currentCoordinatesAgent,
    double rangeDispersal,//buffer around the travel lines forming a polygon where seed could be dispersed
    double mapSize
){//In case you need a reminder of formulas, see bottom of https://fr.wikipedia.org/wiki/Rotation_plane
  double distanceToTarget = distanceTravelledCalculation(
    previousCoordinatesAgent,
    currentCoordinatesAgent
  );
  
  //Insert coordinates as if x axis is travelled path
  // Rcout << distanceToTarget << std::endl;
  NumericVector newTreeLoc = rep(0.0, 2);
  newTreeLoc[0] = runif(2, 0, distanceToTarget)[0];
  // Rcout << newTreeLoc[0] << std::endl;
  newTreeLoc[1] = runif(2, -rangeDispersal, rangeDispersal)[0];
  // Rcout << newTreeLoc[1] << std::endl;
  
  //Backtransform to have them in initial referential system
  NumericVector newTreeLocInitialRef = rep(0.0, 2);
  
  //x coord
  newTreeLocInitialRef[0] = previousCoordinatesAgent[0] + 
    newTreeLoc[0]*(currentCoordinatesAgent[0] - previousCoordinatesAgent[0])/distanceToTarget +
    -newTreeLoc[1]*(currentCoordinatesAgent[1] - previousCoordinatesAgent[1])/distanceToTarget;
    //y coord
    newTreeLocInitialRef[1] = previousCoordinatesAgent[1] + 
    newTreeLoc[0]*(currentCoordinatesAgent[1] - previousCoordinatesAgent[1])/distanceToTarget +
    newTreeLoc[1]*(currentCoordinatesAgent[0] - previousCoordinatesAgent[0])/distanceToTarget;
    
    //Modify while not in the map (because of perceptual range, it can indeed be out of the map)
    while(
      (newTreeLocInitialRef[0] >= mapSize) | (newTreeLocInitialRef[0] <= 0) | //x out of map
        (newTreeLocInitialRef[1] >= mapSize) | (newTreeLocInitialRef[1] <= 0) //y out of map
    ){
      //Rcout << "Out of map" << std::endl;
      double newX = runif(2, 0, distanceToTarget)[0];
      double newY = runif(2, -rangeDispersal, rangeDispersal)[0];
      
      //x coord
      newTreeLocInitialRef[0] = previousCoordinatesAgent[0] + 
        newX*(currentCoordinatesAgent[0] - previousCoordinatesAgent[0])/distanceToTarget +
        -newY*(currentCoordinatesAgent[1] - previousCoordinatesAgent[1])/distanceToTarget;
        //y coord
        newTreeLocInitialRef[1] = previousCoordinatesAgent[1] + 
        newX*(currentCoordinatesAgent[1] - previousCoordinatesAgent[1])/distanceToTarget +
        newY*(currentCoordinatesAgent[0] - previousCoordinatesAgent[0])/distanceToTarget;
    } 
    
    NumericVector testNA_v = newTreeLocInitialRef;
    if(any(is_na(testNA_v))){
      Rcout << "Error coordinates after Dispersal" << std::endl;
      Rcout << currentCoordinatesAgent[1] - previousCoordinatesAgent[1] << " " <<
        currentCoordinatesAgent[0] - previousCoordinatesAgent[0]<< " " <<
          distanceToTarget << std::endl;
    }
    
    //Return the new tree coordinates accounting for Dispersal of seeds. Note, this does not account for the probability of Dispersal, which is to be done indenpently
    return(newTreeLocInitialRef);

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Determine which tree was fruiting
//--------------------

//This function gives a binary vector to tell wether at a given time, a fruit tree is fruiting or not.

// [[Rcpp::export]]
IntegerVector whichInFruit(
    double currentTime,
    NumericVector startTimeVector,//vector of fruiting start date of trees
    NumericVector endTimeVector,//vector of fruiting end date of trees
    double lengthOfCycle
){
  IntegerVector outputVector(startTimeVector.size());
  NumericVector currentTimeVector =  rep(currentTime, startTimeVector.size());
  
  //Translate so 0 matches start Time
  currentTimeVector = currentTimeVector - startTimeVector;
  endTimeVector = endTimeVector - startTimeVector;
  startTimeVector = startTimeVector - startTimeVector;
  
  // //Do the modulo (actually only necessary for timer since there has been a translation)
  currentTimeVector = moduloVector(currentTimeVector, lengthOfCycle);
  endTimeVector = moduloVector(endTimeVector, lengthOfCycle);
  startTimeVector = moduloVector(startTimeVector, lengthOfCycle);
  
  outputVector = ifelse(
    (currentTimeVector >= startTimeVector) & (currentTimeVector <= endTimeVector),
    1,
    0
  );
  
  //Return binary vector that indictes whether tree is in fruit or not 
  return outputVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----- Agent-related

//--------------------
// Move agent
//--------------------

//This function estimates what tree is the best ratio quantity of food/distance based on the agent's knowledge.

// [[Rcpp::export]]
NumericVector moveAgentKnowledgeBased(
    NumericVector currentCoordinatesAgent,
    NumericMatrix treeLoc,
    NumericVector treeFood,//quantity of food at each tree
    std::string whatRule,//Should be "closest" or "farthest"
    int IDcurrentTree
){
  NumericVector outputVector(3);
  //Calculate distance to each tree
  NumericVector vectorDistance = sqrt(pow((treeLoc(_,0) - currentCoordinatesAgent[0]),2) + pow((treeLoc(_,1) - currentCoordinatesAgent[1]),2));
  //Calculate interest for each tree                                      
  NumericVector vectorInterest = floor(treeFood/vectorDistance * 10000000)/10000000;//I used *10000 because with double there are problems with tiny differences that create absence of equality while it should
  
  //Change interest of current location to negative
  
  // //a) determine which tree is current loc: will be the 0 in the vector
  // NumericVector vectorToDetermineCurrentTree = (treeLoc(_,0) - currentCoordinatesAgent[0]) + (treeLoc(_,1) - currentCoordinatesAgent[1]);
  // //b) change interest
  // vectorInterest = ifelse(vectorToDetermineCurrentTree < 0.001, -1000, vectorInterest);
  
  if(IntegerVector::is_na(IDcurrentTree) || IDcurrentTree == -1000){
  }else{
    vectorInterest[IDcurrentTree] = -1000;
  };

  //Select tree(s) with best interest
  IntegerVector whichTreeBest = which_maxRcpp(vectorInterest);
  
  //Rcout << whichTreeBest.size() << std::endl;
  if(whichTreeBest.size() == 0){
    Rcout << "Error, no target despite cognitive foraging" << std::endl;
  }
  if(whichTreeBest.size() > 1){
    NumericVector distanceBestTrees = vectorDistance[whichTreeBest];
    if(whatRule == "closest"){//If several trees, select the closest
      int whichTreeInBestTree(which_minRcpp(distanceBestTrees)[0]);
      int whatTreeIDSelect(whichTreeBest[whichTreeInBestTree]);
      //Rcout << "Closest computed" << std::endl;
      outputVector[0] = treeLoc(whatTreeIDSelect,0);
      outputVector[1] = treeLoc(whatTreeIDSelect,1);
      outputVector[2] = whatTreeIDSelect;
    }else{//If several trees, select the farthest
      int whichTreeInBestTree(which_maxRcpp(distanceBestTrees)[0]);
      int whatTreeIDSelect(whichTreeBest[whichTreeInBestTree]);
      outputVector[0] = treeLoc(whatTreeIDSelect,0);
      outputVector[1] = treeLoc(whatTreeIDSelect,1);
      outputVector[2] = whatTreeIDSelect;
    }
  }else{//If one tree, take it
    outputVector[0] = treeLoc(whichTreeBest[0],0);
    outputVector[1] = treeLoc(whichTreeBest[0],1);
    outputVector[2] = whichTreeBest[0];
  }
  //Return vector of coordinates of the best target (x,y) and the ID of the tree
  return(outputVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Determine whether tree was perceived on the way to the target or not; and at start/end of the movement bout
//--------------------

//This function determines what trees were encountered along the way, considering the agent perceptual range.

// [[Rcpp::export]]
IntegerVector visitedTrees(
    NumericVector previousCoordinatesAgent,
    NumericVector currentCoordinatesAgent,
    NumericMatrix treeLoc,
    double perceptualRange
){
  
  double distanceToTarget = distanceTravelledCalculation(
    previousCoordinatesAgent,
    currentCoordinatesAgent
  );
  
  NumericMatrix newTreeLoc = changeReferentialCoordinatesAlongTravel(
    previousCoordinatesAgent,
    currentCoordinatesAgent,
    treeLoc
  );
  
  NumericVector distanceAtStart = 
    sqrt(
      pow(treeLoc(_,0) - previousCoordinatesAgent[0], 2) +
      pow(treeLoc(_,1) - previousCoordinatesAgent[1], 2) 
    );
    
  NumericVector distanceAtEnd = 
    sqrt(
      pow(treeLoc(_,0) - currentCoordinatesAgent[0], 2) +
      pow(treeLoc(_,1) - currentCoordinatesAgent[1], 2) 
    );
    
  //Determine trees that were visited en route, at start or end of movement bout
  IntegerVector isSeenVector = ifelse(
    (
      (
      (abs(newTreeLoc(_,1)) <= perceptualRange) &//The tree encountered en route is no more than visionRadius away from the linear path to the target
      (newTreeLoc(_,0) >= 0) & //This tree is in the same direction than the target
      (newTreeLoc(_,0) <= distanceToTarget)
      ) |// This tree is not farther than the target
      (distanceAtStart <= perceptualRange) | //Tree was perceived at start
      (distanceAtEnd <= perceptualRange)
    ), //Tree was perceived at end
      1,
      0
  );
  
  //Return a binary vector to indicated which trees were perceived during travel
  return(isSeenVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Determine whether tree was perceived on the way to the target or not; and at start/end of the movement bout
//--------------------

// [[Rcpp::export]]
NumericVector Rcpp_sort(NumericVector x, NumericVector y) {
  // Order the elements of x by sorting y
  // First create a vector of indices
  IntegerVector idx = seq_along(x) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  // And return x in that order
  return(x[idx]);
}

// [[Rcpp::export]]
IntegerVector visitedTreesSequence(NumericMatrix treeLoc, 
                                   IntegerVector isSeenVector,//binary vector of trees if visited (1) or not (0)
                                   NumericVector currentLoc){
  //Distance of from current location to every visited tree
  NumericVector distanceVector =
    pow((treeLoc(_,0) - currentLoc[0]), 2) +
    pow((treeLoc(_,1) - currentLoc[1]), 2);
  //Create matrix of tree with distance and id
  IntegerVector idVector = seq_along(treeLoc(_,0)) - 1;
  NumericMatrix distanceMatrix(idVector.size(), 2); 
  distanceMatrix(_, 0) = distanceVector;
  distanceMatrix(_, 1) = idVector;
  //Subset only matrix of distances of visited trees
  NumericMatrix distanceMatrix_rdc(sum(isSeenVector), 2);
 
  int count = 0;
  for(int i = 0; i < treeLoc(_,0).size(); i++){
    if(isSeenVector[i] == 1){
      distanceMatrix_rdc(count, 0) = distanceMatrix(i, 0);
      distanceMatrix_rdc(count, 1) = distanceMatrix(i, 1);
      count += 1;
      // Rcout << count << std::endl;
    }
  }
  //Reorder in function of distance to get sequence of visit
  NumericVector outputVector = Rcpp_sort(distanceMatrix_rdc(_, 1), distanceMatrix_rdc(_, 0));
  return(IntegerVector(outputVector));
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------


//--------------------
// Determine which tree is memorised
//--------------------

//This functions outputs a matrix of ones and zeros in two columns to tell whether spatial (first col) or temporal (second col) is known.

// [[Rcpp::export]]
IntegerMatrix whichDistanceAndFoodKnown(
  int numberTrees,
  double spatialKnowledgeRate,
  double temporalKnowledgeRate
){
  IntegerMatrix outputMatrix(numberTrees, 2);
  IntegerVector IDTreeVector = seq_len(numberTrees) - 1;

  //Check that temporal knowledge is not more than spatial knowledge
  if(spatialKnowledgeRate < temporalKnowledgeRate){
    Rcpp::stop("Problem, spatial knowledge is less than temporal knowledge");
    return outputMatrix;
  } else{

    // First sample trees for which space is known
    if(floor(spatialKnowledgeRate * numberTrees) > 0.01){
      IntegerVector spaceKnownID = sample(IDTreeVector, floor(spatialKnowledgeRate * numberTrees));
      for(int i = 0; i < spaceKnownID.size(); i++){
        outputMatrix(spaceKnownID[i], 0) = 1;
        //Rcout << IDTreeVector[i] << std::endl;
        //Rcout << spaceKnownID[i] << std::endl;
      }
      //Second, sample trees for which time is known. Are included in those for which space is known.
      if(floor(temporalKnowledgeRate/spatialKnowledgeRate * spaceKnownID.size()) > 0.01){
        IntegerVector timeKnownID = sample(spaceKnownID, floor(temporalKnowledgeRate/spatialKnowledgeRate * spaceKnownID.size()));
        for(int i = 0; i < timeKnownID.size(); i++){
          outputMatrix(timeKnownID[i], 1) = 1;
        }
      }
    }
    //Return binary matrix with 1 indicating known trees from a spatial point of view (first column), or temporal point of view (column 2) 
    return outputMatrix;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Determine which tree is perceived
//--------------------

//This functions estimates what are the tree in the agent's perceptual range based on its location.

// [[Rcpp::export]]
IntegerVector whichPerceived(
  NumericMatrix treeLoc,
  NumericVector currentLoc,
  double perceptualRange
){
  IntegerVector outputVector(treeLoc(_,0).size());
  
  NumericVector distanceToTree = 
    sqrt(
      pow(treeLoc(_,0) - currentLoc[0], 2) +
      pow(treeLoc(_,1) - currentLoc[1], 2)
    );
  
  outputVector = ifelse(
    distanceToTree <= perceptualRange,
    1,
    0
  );
  
  //Return binary vector to indicate which trees are perceived from position
  return outputVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Randomwalk
//--------------------

// [[Rcpp::export]]
NumericVector moveAgentRandom(
  NumericVector currentLocation,
  double mapSize,
  double exponentialRate
){
  NumericVector outputVector(2);
  
  double stepLength = ceil(rexp(1, exponentialRate)[0]);
  double headingInRadian = runif(1,-atan(1)*4, atan(1)*4)[0];//Reminder: atan(1)*4 == pi 
  double newX(stepLength * cos(headingInRadian) + currentLocation[0]);
  double newY(stepLength * sin(headingInRadian) + currentLocation[1]);
  //Rcout << newX << " " << newY << std::endl;
  
  while((newX < 0) | (newX > mapSize) | (newY < 0) | (newY > mapSize)){
    stepLength = ceil(rexp(1, exponentialRate)[0]);
    headingInRadian = runif(1, -atan(1)*4, atan(1)*4)[0];//Reminder: atan(1)*4 == pi 
    newX = cos(headingInRadian) * stepLength + currentLocation[0];
    newY = sin(headingInRadian) * stepLength + currentLocation[1];
    //Rcout << newX << " " << newY << std::endl;
  }
  
  outputVector[0] = newX;
  outputVector[1] = newY;

  //Return vector of coordinates
  return outputVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Filter unknown trees
//--------------------

//This function transform the matrix of tree location + phenology based on the agent's knowledge.
//In particular, distance of unknown trees is faked to be extremly high, so they will never be chosen.
//For temporal knowledge, a fake value is added.

// [[Rcpp::export]]
NumericMatrix filterKnowledge(
  NumericMatrix treeLoc,
  NumericVector foodTree,
  IntegerMatrix whichDistanceAndFoodKnown,
  NumericVector timeLastVisit_v,
  double currentTime,
  double noReturnTime,
  double whatValueUnknownTemporal,
  double mapSize
){
  NumericMatrix outputMatrix(timeLastVisit_v.size(), 3);
  outputMatrix(_, 0) = treeLoc(_,0);
  outputMatrix(_, 1) = treeLoc(_,1);
  outputMatrix(_, 2) = foodTree;
  
  for(int i = 0; i < timeLastVisit_v.size(); i++){
    
    //Filtrate coordinates if unknown spatially
    if(whichDistanceAndFoodKnown(i,0) == 0){
      outputMatrix(i, 0) = mapSize * 10000000;
      outputMatrix(i, 1) = mapSize * 10000000;
      //Filtrate coordinates if unknown temporally
      if(whichDistanceAndFoodKnown(i,1) == 0){
        outputMatrix(i, 2) = -1000;//If no temporal knowledge, trees will be used only if no other trees with temporal knowledge yield fruit
      }
    }else{
      //Filtrate coordinates if unknown temporally
      if(whichDistanceAndFoodKnown(i,1) == 0){
        outputMatrix(i, 2) = whatValueUnknownTemporal;//If no temporal knowledge, trees will be used only if no other trees with temporal knowledge yield fruit
      }
    }
    
    //Filtrate among trees recently visited
    if(timeLastVisit_v[i] >= currentTime - noReturnTime){
       //& 
       // whichDistanceAndFoodKnown(i,0) == 1 &
       // whichDistanceAndFoodKnown(i,1) == 0
      outputMatrix(i, 2) = -1000;
    }else{}
  }

  //Return a matrix with new coordinates and food to use (i.e. unknown replaced by NAs)
  return outputMatrix;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Modify last visit time, last consumption time, cancel food consumed at end of production cycle
//--------------------

//This function returns a matrix of last visit, last visit with food, tot food eaten for a given year (i.e. will reset info each end of cycle for the tree, based on end fruiting date).
//This is meant for each move to reupdate this info.

// [[Rcpp::export]]
NumericMatrix updateVisitTimesAndConsumptionCycle(
    double currentTime,
    NumericVector startTimeVector,
    NumericVector endTimeVector,
    double lengthOfCycle,
    NumericVector foodConsumedVector,
    IntegerVector visitedTrees_v,
    NumericVector lastTimeVisitTrees_v,
    NumericVector lastTimeVisitTreesWithFruit_v,
    NumericVector foodQuantityAtTree_v
){
  NumericMatrix outputMatrix(startTimeVector.size(), 3);
  NumericVector currentTimeVector =  rep(fmod(currentTime + 10*lengthOfCycle, lengthOfCycle), startTimeVector.size());
  
  //Translate so 0 matches start Time
  currentTimeVector = currentTimeVector - startTimeVector;
  endTimeVector = endTimeVector - startTimeVector;
  startTimeVector = startTimeVector - startTimeVector;
  
  // //Do the modulo (actually only necessary for timer since there has been a translation)
  currentTimeVector = moduloVector(currentTimeVector, lengthOfCycle);
  endTimeVector = moduloVector(endTimeVector, lengthOfCycle);
  startTimeVector = moduloVector(startTimeVector, lengthOfCycle);
  
  //Modify last time visit
  outputMatrix(_, 0) = ifelse(visitedTrees_v == 1, currentTime, lastTimeVisitTrees_v);
  
  //Modify last time feeding
  outputMatrix(_, 1) = ifelse((visitedTrees_v == 1) & (currentTimeVector > startTimeVector) & (currentTimeVector < endTimeVector),
         currentTime,
         lastTimeVisitTreesWithFruit_v);
  
  //Modify food consumed
  //a) due to recent movement
  outputMatrix(_, 2) = ifelse(
         visitedTrees_v == 1,
         foodConsumedVector + foodQuantityAtTree_v,
         foodConsumedVector);
  
  //b) if new season for tree begins, i.e. put 0 if fruiting period ends
  outputMatrix(_, 2) = ifelse(currentTimeVector > endTimeVector,
         0,
         outputMatrix(_, 2));
  
  //Return an output matrix: Col 1 = Time last visit, Col 2 = Time last visit eating, Col 3 = foodConsumed at tree
  return outputMatrix;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Nearest neighbour mean and sd
//--------------------

//This function calculates the distance to all trees.

// [[Rcpp::export]]
NumericVector nearestNeighbourDistance(
  NumericMatrix locMatrix,
  int maxDistancePossible
){
  int numberObs = locMatrix(_,1).size();
  NumericVector distanceToNN = rep(0.0, numberObs);
  for(int loop = 0; loop < numberObs; loop++){
    //Distance to all locs
    NumericVector distanceToLocs = pow(locMatrix(_, 0) - locMatrix(loop, 0), 2) + pow(locMatrix(_, 1) - locMatrix(loop, 1), 2);
      
    //Consider itself as max possible distance, to avoid counting it as closest target
    distanceToLocs[loop] = maxDistancePossible;
    
    //Save minimal distance
    distanceToNN[loop] = min(distanceToLocs);
  }

  //Return a vector with the mean and sd of distance to nearestNeighbour
  NumericVector outputVector = NumericVector::create(mean(distanceToNN), sd(distanceToNN));
  return outputVector;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------
// Lloyd index of patchiness
//--------------------

//This function calcules the patchiness index developped by Lloyd.

// [[Rcpp::export]]
NumericVector lloydIndex(
    NumericMatrix treeLoc,
    int mapSize,
    int quadratSize
){
  NumericMatrix matrixCount(mapSize/quadratSize, mapSize/quadratSize);
  //Count trees in each quadrat
  for(int tree = 0; tree < treeLoc(_,1).size(); tree++){
    int x = floor(treeLoc(tree,0)/quadratSize);
    int y = floor(treeLoc(tree,1)/quadratSize);
    matrixCount(x,y) += 1;
    //Rcout << matrixCount(x,y) << " " << x << " " << y << std::endl;
  }
  //Calculate indices
  
  //a) Transform matrix in a vector
  IntegerVector matrixTransformedToVector(mapSize/quadratSize*mapSize/quadratSize);
  for(int counter = 0; counter < mapSize/quadratSize; counter++){
    IntegerVector toModifyVector = seq(0, mapSize/quadratSize - 1) + mapSize/quadratSize*counter;
    for(int counter2 = 0; counter2 < toModifyVector.size(); counter2++){
      //Rcout << matrixCount(counter2,counter) << std::endl;
      matrixTransformedToVector[toModifyVector[counter2]] = matrixCount(counter2,counter);
    }
  }
  //b) compute calculations
  NumericVector output_v(2);
  double meanValueVector = mean(matrixTransformedToVector);
  IntegerVector matrixTransformedToVectorMinus1 = matrixTransformedToVector - 1;
  IntegerVector productVectorVectorMinus1 = matrixTransformedToVector * matrixTransformedToVectorMinus1;

  double numerator = sum(productVectorVectorMinus1);
  double denominator = sum(matrixTransformedToVector);
  //Rcout << numerator << " " << denominator << std::endl;
  output_v[0] = numerator/denominator;
  output_v[1] = output_v[0]/meanValueVector;
  
  //Return mean crowding and index of patchiness
  return output_v;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]

// This function will return a boolean vector telling wether the vectors have same values or not. Vectors should be of same size.

bool compareVectors(NumericVector Vector1, NumericVector Vector2){
  bool output = any(ifelse(Vector1 == Vector2, true, false)).is_true();
  return output;
}

//--------------------
// Run simulation
//--------------------

//This function runs the simulation in which a forager with potentially knowledge forage in a closed environment alone.

// [[Rcpp::export]]
void runSimulation(
const int cycleLimitNumber, 
const int repetitionNumber, 
const double timeDelayForDispersal,
const double torporTime, 
const bool saveTreeMap, 
IntegerVector samplingMapTime_v,
std::string nameInit,
const int mapSize, 
const int quadratSize,
const int numberTrees,
const bool homogeneousDistribution,
const int treeClusterNumber,
const int treeClusterSpread,
NumericMatrix treeLocInit_m,
NumericMatrix fruitingTimesInit_m,
const NumericVector maximumFoodToYield_v, 
const double cycleLength, 
const int fruitingLength,
const double noReturnTime, 
const double whatValueUnknownTemporal, 
std::string whatRule, 
const double exponentialRate,
double perceptualRange,
double spatialKnowledgeRate,
double temporalKnowledgeRate,
const double speed, 
const double DispersalProbability,
const bool useProvidedMap,
double intensityCompetitionForSpace = 0.85,
bool moveOnlyToFruitingTrees = false,
bool moveOnlyToTarget = false,
bool linear = true,
bool learning = false
){
  ////----------------------------------------------------
  //// INITIALISATION
  ////----------------------------------------------------

  set_seed(42);
  // for(int loopingForMapSave = 0; loopingForMapSave < numberTrees; loopingForMapSave++){
  //   Rcout << loopingForMapSave << " " << treeLocInit_m(loopingForMapSave, 0) << " " <<  treeLocInit_m(loopingForMapSave, 1) << std::endl;
  // }
  
  //---- General
  //List of parameters

  double timer(0.0);
  double timerPrevious(timer);
  int eventsWithNoFood(0);

  //File for final record of foraging efficiency (cumulated)
  std::string pathOfTheFileOutputInit(nameInit);
  pathOfTheFileOutputInit.append("_p");
  pathOfTheFileOutputInit.append(std::to_string(perceptualRange));
  pathOfTheFileOutputInit.append("_s");
  pathOfTheFileOutputInit.append(std::to_string(spatialKnowledgeRate));
  pathOfTheFileOutputInit.append("_t");
  pathOfTheFileOutputInit.append(std::to_string(temporalKnowledgeRate));
  pathOfTheFileOutputInit.append("_r");
  pathOfTheFileOutputInit.append(std::to_string(repetitionNumber));

  //File for routine records (i.e. sequence of trees)
  std::string pathOfTheFileOutputRoutine = pathOfTheFileOutputInit;
  pathOfTheFileOutputRoutine.append("_Routine.txt");
  std::ofstream outputRoutine;
  //Open connection to output file
  outputRoutine.open(pathOfTheFileOutputRoutine.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
  
  //File for output of continuous foraging efficiency monitoring
  std::string pathOfTheFileOutputInitEfficiencyContinuous = pathOfTheFileOutputInit;
  pathOfTheFileOutputInitEfficiencyContinuous.append("_EfficiencyContinuous.txt");
  std::ofstream outputFluxEfficiencyContinuous;
  //Open connection to output file
  outputFluxEfficiencyContinuous.open(pathOfTheFileOutputInitEfficiencyContinuous.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
  if(outputFluxEfficiencyContinuous)  // Write the file only if correctly opened
  {
    outputFluxEfficiencyContinuous <<
      "Time" << " " <<
        "Efficiency" << std::endl;//" " <<
          // "Nearest_neighbour_mean" << " " <<
          //   "Nearest_neighbour_sd" << " " <<
          //     "meanCrowding" << " " <<
          //       "patchinessIndex" << std::endl;
  }else //If impossible to write in the file, say it
  {
    Rcout << "ERROR: Impossible to open the output file about continuous efficiency" << std::endl;
  }

  int hasSavedFinalResults(0);
  
  //---- Environment related
  int lowerBorder(0); //Lower border of a coordinate (x or y)
  int upperBorder(mapSize); //Upper border of a coordinate (x or y)
  
  NumericMatrix treeLocReal_m = treeLocInit_m;
  
  double spaceTreeSq = intensityCompetitionForSpace * intensityCompetitionForSpace * mapSize * mapSize / (numberTrees * 3.14159265358979323846);//(perceptualRange - 0.10 * perceptualRange)*(perceptualRange - 0.10 * perceptualRange);//(perceptualRange/2)*(perceptualRange/2);//The space expected to be occupied by the tree: considering a circular buffer, area of non-overlapping trees should be equal to that of the map//
  //Rcout << spaceTreeSq << std::endl;
  //Change according to whether to use provided map or not
  if(useProvidedMap){
    //Let initial settings;
  }else{
      //Distribute trees
      treeLocInit_m = distributionTree(
        numberTrees,
        lowerBorder,
        upperBorder,
        true,
        treeClusterNumber,
        treeClusterSpread
      );
      //Resample fruiting dates
      fruitingTimesInit_m = dateTree(
        numberTrees,
        cycleLength,
        fruitingLength
      ); 
    }

  NumericMatrix fruitingTimes_m = fruitingTimesInit_m;
  
  //Initialise current food available
  NumericVector foodConsumed_v = rep(0.0, numberTrees);//Null vector at first
  //NumericVector foodConsumedPrevious_v = rep(0.0, numberTrees);//Null vector at first
  NumericVector foodQuantityAtTree_v = foodTree(
    0,
    fruitingTimesInit_m(_,0),
    fruitingTimesInit_m(_,1),
    cycleLength,
    foodConsumed_v,
    maximumFoodToYield_v,
    linear
  );
  
  //---- Agent-related
  //Initialise location
  NumericVector agentCurrentLocation = NumericVector::create(runif(1, 0, mapSize)[0], runif(1, 0, mapSize)[0]);
  //Rcout << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
  NumericVector agentPreviousLocation = agentCurrentLocation;

  //Initialise visited trees
  NumericVector lastTimeVisitTrees_v = rep(0.0, numberTrees);
  NumericVector lastTimeVisitTreesWithFruit_v = rep(0.0, numberTrees);
  IntegerVector visitedTrees_v = rep(0, numberTrees);//Null vector
  IntegerVector disseminated_v = rep(1, numberTrees);//Null vector
  NumericVector lastTimeDispersal_v = rep(0.0, numberTrees);
  
  //Initialise knowledge on trees
  IntegerMatrix spatiotemporalKnowledge_m = whichDistanceAndFoodKnown(
    numberTrees,
    spatialKnowledgeRate,
    temporalKnowledgeRate
  );

  NumericMatrix treeLoc_knowledge_m = filterKnowledge(
    treeLocReal_m,
    foodQuantityAtTree_v,
    spatiotemporalKnowledge_m,
    lastTimeVisitTrees_v,
    timer,
    noReturnTime,
    whatValueUnknownTemporal,
    mapSize
  );

  // Rcout << "Food knowledge: " << sum(ifelse(treeLoc_knowledge_m(_,2) >= 0, treeLoc_knowledge_m(_,2), 0)) << std::endl;
  // 
  //Results: distance travelled and food consumed,
  double totDistanceTravelled(0);
  double totFoodEaten(0);
  int samplingMapCounter(0);
  int IDcurrentTree = NA_INTEGER;
  int IDpreviousTree = NA_INTEGER;
  //bool compVec = compareVectors(treeLocInit_m(_,1), treeLocReal_m(_,1));
  
  ////----------------------------------------------------
  //// RUN SIMULATION
  ////----------------------------------------------------
  
  while(timer < (cycleLimitNumber * cycleLength + 5 * cycleLength)){//last five cycles will be to analyse routine
    //Rcout << timer << std::endl;
    
    //Checking fruiting dates
    // if(timer >= 7206){
    //   Rcout << fruitingTimes_m(_,0).size() << std::endl;
    //   for(int mat = 0; mat < fruitingTimes_m(_,0).size(); mat++){
    //     if(abs(fruitingTimes_m(mat,1) - fruitingTimes_m(mat,0)) > fruitingLength + 1 || abs(fruitingTimes_m(mat,1) - fruitingTimes_m(mat,0)) < fruitingLength - 1){
    //       Rcout << "ERROR WITH FRUITING DATES TREE " << mat << std::endl;
    //     }
    //   }
    // }
    
    // compVec = compareVectors(treeLocInit_m(_,1), treeLocReal_m(_,1));
    // if(compVec){
    //   Rcout << "Comparison tree loc with init 1: DIFFERENT " << timer << std::endl;
    // }else{
    //   Rcout << "Comparison tree loc with init 1: SAME " << timer << std::endl;
    // }
    // 
    
    //Update knowledge matrix
    //a) Discard depletion, which is unknown
    // NumericVector nullVector = rep(0.0, numberTrees);
    // NumericVector foodQuantityAtTree_nodepletion_v =
    //   foodTree(
    //     timer,
    //     fruitingTimes_m(_,0),
    //     fruitingTimes_m(_,1),
    //     cycleLength,
    //     nullVector,
    //     maximumFoodToYield_v,
    //     linear
    //   );

    //b) Update matrix: use foodQuantityTree_nodepletion_v if you don't want to account for depletion
    treeLoc_knowledge_m = filterKnowledge(
      treeLocReal_m,
      foodQuantityAtTree_v,
      spatiotemporalKnowledge_m,
      lastTimeVisitTrees_v,
      timer,
      noReturnTime,
      whatValueUnknownTemporal,
      mapSize
    );

    //Adjust knowledge based on perception
    
    IntegerVector whichTreeInPerceptualRange = whichPerceived(
      treeLocReal_m,
      agentCurrentLocation,
      perceptualRange
    );
    
    //Rcout << "Knowledge spatial and temporal " << sum(spatiotemporalKnowledge_m(_,0)) << " " << sum(spatiotemporalKnowledge_m(_,1)) << std::endl;
    //Rcout << "Food in vision " << sum(ifelse(whichTreeInPerceptualRange == 1, foodQuantityAtTree_v, 0)) << std::endl; //Food in perception: should be 0 after first move
    //Rcout << "Food knowledge prior vision " << sum(treeLoc_knowledge_m(_,2)) << std::endl; //Food in perception: should be 0 after first move
    
    for(int mat = 0; mat < whichTreeInPerceptualRange.size(); mat++){
      if((whichTreeInPerceptualRange[mat] == 1) && 
         (foodQuantityAtTree_v[mat] > 0.01) &&
         (lastTimeVisitTrees_v[mat] <= (timer - noReturnTime))
        ){
          treeLoc_knowledge_m(mat,0) = treeLocReal_m(mat,0);
          treeLoc_knowledge_m(mat,1) = treeLocReal_m(mat,1);
          treeLoc_knowledge_m(mat,2) = foodQuantityAtTree_v[mat];//I use > 0.01 to not change those recently visited
      }
    }

    NumericVector targetCoordinates = NumericVector::create(0,0);//runif(2, 0, mapSize);

    //*****************
    // 1) MOVE
    //*****************

    int targetIDInit = -1000;

    //Enter torpor if no food, otherwise move
    bool enterTorpor(false);
    double distanceTravelled(0);
    if(sum(foodQuantityAtTree_v) <= 1){
      enterTorpor = true;
    }

    if(enterTorpor){
      //Rcout << "Torpor" << std::endl;
      timer += torporTime;
      eventsWithNoFood += 1;
      //Rcout << timer << std::endl;
    }else{
      //Move to target
      //a) Chose target
      agentPreviousLocation = agentCurrentLocation;

      //Calculate amount of food estimated: depletion not accounted because recently visited trees are not visited again, otherwise they are considered as unvisited, thus with no depletion

      //c)Move agent randomly if no food known/perceived
      // NumericVector targetCoordinates = agentCurrentLocation;
      // if(timer >= 7206){
      //   Rcout << "Food quantity " << sum(foodQuantityAtTree_v) << std::endl;
      //   Rcout << "Food knowledge " << sum(ifelse(treeLoc_knowledge_m(_,2) > 0.001, treeLoc_knowledge_m(_,2), 0)) << std::endl;
      //   Rcout << "Food knowledge " << sum(ifelse(treeLoc_knowledge_m(_,2) > 0.001, 1, 0)) << std::endl;
      //   Rcout << "ID max food " << which_maxRcpp(treeLoc_knowledge_m(_,2))[0] << std::endl;
      // }
      bool forageBasedOnMemory(false);
      if(sum(ifelse(treeLoc_knowledge_m(_,2) > 0.001, 1, 0)) > 0){
        forageBasedOnMemory = true;
      }

      // if(timer >= 7206){
      //   Rcout << forageBasedOnMemory << std::endl;
      //   Rcout << "Max food is " << max(treeLoc_knowledge_m(_,2)) << std::endl;
      //   Rcout << "Max food tree id " << which_maxRcpp(treeLoc_knowledge_m(_,2)) << std::endl;
      // }

      if(forageBasedOnMemory){
        // Rcout << "I move based on cognition" << std::endl;

        //Select the initial long distance target
        NumericVector outputMovement = moveAgentKnowledgeBased(
         agentCurrentLocation,
         treeLoc_knowledge_m(_,Range(0,1)),
         treeLoc_knowledge_m(_,2),
         "closest",//Should be "closest" or "farthest"
         IDcurrentTree = IDcurrentTree
       );
        
        targetCoordinates[0] = outputMovement[0];
        targetCoordinates[1] = outputMovement[1];
        targetIDInit = outputMovement[2];

        // if(timer >= 7206){
        //   Rcout << "Time is " << timer << " Last time visit is " << lastTimeVisitTrees_v[outputMovement[2]] << std::endl;
        //   Rcout << "Current tree is " << IDcurrentTree << std::endl;
        //   Rcout << "Target is tree " << outputMovement[2] << std::endl;
        //   Rcout << "Food at target is " << foodQuantityAtTree_v[outputMovement[2]] << std::endl;
        //   Rcout << "No return time " << noReturnTime << " Target last visited " << timer - lastTimeVisitTrees_v[outputMovement[2]] << std::endl;
        // }

        if(moveOnlyToTarget){
          //Rcout << "Moving target" << std::endl;
          IDpreviousTree = IDcurrentTree;
          IDcurrentTree = targetIDInit;
        }else{
          //Rcout << "Moving not target" << std::endl;
          //Assess all the trees on the way
          IntegerVector expectedVisitedTrees_v = visitedTrees(
            agentCurrentLocation,
            targetCoordinates,
            treeLocReal_m,
            perceptualRange
          );
          // if(timer >= 7206){
          //   Rcout << "Worked up here" << std::endl;
          // }

          //Assess all the trees which were visited sufficiently a long time ago
          IntegerVector visitedLongTimeAgo_v = rep(0, expectedVisitedTrees_v.size());
          for(int v = 0; v < visitedLongTimeAgo_v.size(); v++){
            if(lastTimeVisitTrees_v[v] <= (timer - noReturnTime)){
              visitedLongTimeAgo_v[v] = 1;
            }
          }
          // if(timer >= 7206){
          //   Rcout << "Worked up here 2" << std::endl;
          // }
          //Reconsider trees because of no return memory
          expectedVisitedTrees_v = expectedVisitedTrees_v * visitedLongTimeAgo_v;

          // if(timer >= 7206){
          //   Rcout << "Worked up here 3" << std::endl;
          // }

          //Assess all the trees with food
          IntegerVector treesWithFood_v = rep(0, expectedVisitedTrees_v.size());
          if(moveOnlyToFruitingTrees){
            for(int v = 0; v < treesWithFood_v.size(); v++){
              if(foodQuantityAtTree_v[v] > 0.001){
                treesWithFood_v[v] = 1;
              }
            }
            //Reconsider trees as will be visited if had food
            expectedVisitedTrees_v = expectedVisitedTrees_v * treesWithFood_v;
          }
          //Have the theoretical sequence of visits
          IntegerVector whatVisitedTrees_v = visitedTreesSequence(
            treeLocReal_m,
            expectedVisitedTrees_v,
            agentCurrentLocation
          );
          //if(timer >= 364.8 & timer < 365.3){
          //   // Rcout << "Target status: seen  " << expectedVisitedTrees_v[outputMovement[2]] << std::endl;
          //   // Rcout << "Target status: visited long time " << visitedLongTimeAgo_v[outputMovement[2]] << std::endl;
          //   // Rcout << "Target status: with food " << treesWithFood_v[outputMovement[2]] << std::endl;
          // Rcout << "Current position is " << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
          //   for(int visit = 0; visit < whatVisitedTrees_v.size(); visit++){
          //     Rcout << "Tree to be visited in order: " << whatVisitedTrees_v[visit] << std::endl;
          //   }
          //    Rcout << "Initial target was " << outputMovement[2] << std::endl;
          //    Rcout << "Last tree expected to be visited " << whatVisitedTrees_v[whatVisitedTrees_v.size() - 1] << std::endl;
          //    Rcout << "First tree expected to be visited " << whatVisitedTrees_v[0] << std::endl;
          //    Rcout << "Was visited long time ago " << visitedLongTimeAgo_v[whatVisitedTrees_v[0]] << std::endl;
          //    Rcout << "Time ast visit " << lastTimeVisitTrees_v[whatVisitedTrees_v[0]] << std::endl;
          //    Rcout << "Size of vector visited" << whatVisitedTrees_v.size() << std::endl;
          //}

          // if(timer >=2 && timer <= 2.1){
          //   Rcout << "Agent loc " << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
          //   Rcout << "306 " << treeLocReal_m(306,0) << " " << treeLocReal_m(306,1) << " " << foodQuantityAtTree_v [306] << std::endl;
          //   Rcout << "450 " << treeLocReal_m(450,0) << " " << treeLocReal_m(450,1) << " " << foodQuantityAtTree_v [450] << std::endl;
          // }
               
          if(whatVisitedTrees_v.size() > 0){
            // Consider target as the closest tree with food that should have been seen en route to the target
            IDpreviousTree = IDcurrentTree;
            IDcurrentTree = whatVisitedTrees_v[0];
            targetCoordinates[0] = treeLocReal_m(whatVisitedTrees_v[0],0);
            targetCoordinates[1] = treeLocReal_m(whatVisitedTrees_v[0],1);
          }else{
            //Here, the target can eventually be less than 0.001, as it is the best ratio which is chosen, so the "else" takes care of it, as the target is now "not expected to be visited"
            IDpreviousTree = IDcurrentTree;
            IDcurrentTree = outputMovement[2];
            targetCoordinates[0] = treeLocReal_m(outputMovement[2],0);
            targetCoordinates[1] = treeLocReal_m(outputMovement[2],1);
          }
        }
      }else{
        //Rcout << "I move randomly" << std::endl;
        targetCoordinates = moveAgentRandom(
          agentCurrentLocation,
          mapSize,
          exponentialRate
        );
        if(moveOnlyToTarget){
          IDpreviousTree = IDcurrentTree;
          IDcurrentTree = -1000;//NA_INTEGER;
        }else{
          IntegerVector expectedVisitedTrees_v = visitedTrees(
            agentCurrentLocation,
            targetCoordinates,
            treeLocReal_m,
            perceptualRange
          );
          //Assess all the trees which were visited sufficiently a long time ago
          IntegerVector visitedLongTimeAgo_v = rep(0, expectedVisitedTrees_v.size());
          for(int v = 0; v < visitedLongTimeAgo_v.size(); v++){
            if(lastTimeVisitTrees_v[v] <= (timer - noReturnTime)){
              visitedLongTimeAgo_v[v] = 1;
            }
          }

          expectedVisitedTrees_v = expectedVisitedTrees_v * visitedLongTimeAgo_v;

          //Assess all the trees with food
          IntegerVector treesWithFood_v = rep(0, expectedVisitedTrees_v.size());
          if(moveOnlyToFruitingTrees){
            for(int v = 0; v < treesWithFood_v.size(); v++){
              if(foodQuantityAtTree_v[v] > 0.001){
                treesWithFood_v[v] = 1;
              }
            }
            //Reconsider trees as will be visited if had food
            expectedVisitedTrees_v = expectedVisitedTrees_v * treesWithFood_v;
          }
          //Have the theoretical sequence of visits
          IntegerVector whatVisitedTrees_v = visitedTreesSequence(
            treeLocReal_m,
            expectedVisitedTrees_v,
            agentCurrentLocation
          );
          //If trees of interest along path, stop there
          if(whatVisitedTrees_v.size() > 0){
            IDpreviousTree = IDcurrentTree;
            IDcurrentTree = whatVisitedTrees_v[0];
            targetCoordinates[0] = treeLocReal_m(whatVisitedTrees_v[0],0);
            targetCoordinates[1] = treeLocReal_m(whatVisitedTrees_v[0],1);
          }else{
            IDpreviousTree = IDcurrentTree;
            IDcurrentTree = -1000;//NA_INTEGER;
          }
        }
      }

      // if(timer >=2 && timer <= 10 && IDcurrentTree  != -1000){
      //   Rcout << timer << " ID Target " << IDcurrentTree << std::endl;
      //   Rcout << "Food previous " << foodQuantityAtTree_v[IDpreviousTree] << std::endl;
      //   Rcout << "Food target " << foodQuantityAtTree_v[IDcurrentTree] << std::endl;
      //   Rcout << "Food knowledge target " << treeLoc_knowledge_m(IDcurrentTree,2) << std::endl;
      //   Rcout << "Food knowledge target " << max(treeLoc_knowledge_m(_,2)) << std::endl;
      //   Rcout << "Coordinates agent " << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
      //   Rcout << "Coordinates target " << targetCoordinates[0] << " " << targetCoordinates[1] << std::endl;
      // }

      // UPDATE ROUTINE
      if(timer >= cycleLimitNumber * cycleLength){
        if(enterTorpor){

        }else{
          //Save routine: it will be when dispersal does not occur
          if(outputRoutine){
            // Write the file only if correctly opened
            //Add some labels to be able to retrieve: points revisited because random movements
            if(forageBasedOnMemory){
              if(IDcurrentTree == targetIDInit){
                outputRoutine << "T" << IDcurrentTree << " timer " << timer << std::endl;
              }else{
                outputRoutine << "t" << IDcurrentTree << " timer " << timer << std::endl;//Now moving from tree to tree progressively, should not be present anymore
              }
            }else{
                outputRoutine << "r" << IDcurrentTree << " timer " << timer << std::endl;
              }
            }
          }
          // if(timer >= 7206){
          //   Rcout << "Routine updated" << std::endl;
          // }
        }

      agentCurrentLocation = targetCoordinates;
      // if(timer >= 365 & timer < 365.1){
      //   Rcout << "Tree target " << IDcurrentTree << std::endl;
      //   Rcout << "Previous tree " << IDpreviousTree << std::endl;
      //   Rcout << "Location target " << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
      // }

      //b) Update visited trees in memory and food depletion

      // The following rule was only valid in case of telescopic arms
      // visitedTrees_v = visitedTrees(
      //   agentPreviousLocation,
      //   agentCurrentLocation,
      //   treeLocReal_m,
      //   perceptualRange
      // );
      // Now this is the true rule: is moving from tree to tree (either just to the target, all trees encountered on the way, just the fruiting trees encountered on the way)
      visitedTrees_v = rep(0, visitedTrees_v.size());
      if(IntegerVector::is_na(IDcurrentTree) || IDcurrentTree == -1000){
        
      }else{
        visitedTrees_v[IDcurrentTree] = 1;
      }
      if(IntegerVector::is_na(IDpreviousTree) || IDpreviousTree == -1000){
        
      }else{
        visitedTrees_v[IDpreviousTree] = 1;
      }

      // Rcout << sum(visitedTrees_v) << " trees were visited at this round!" << std::endl;
      // Rcout << visitedTrees_v[IDcurrentTree] << std::endl;

      //Update Distance travelled and timing
      distanceTravelled = distanceTravelledCalculation(
        agentPreviousLocation,
        agentCurrentLocation
      );
      
      // if(timer >= 7206){
      //   Rcout << "Distance travelled of " << distanceTravelled  << std::endl;
      // }
      if(distanceTravelled <= 0){
        
        std::string pathOfTheFileOutputError = pathOfTheFileOutputInit;
        pathOfTheFileOutputError.append("_Error.txt");
        std::ofstream outputFluxError;
        //Open connection to output file
        outputFluxError.open(pathOfTheFileOutputError.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
        outputFluxError << "Time is " << timer << std::endl;
        outputFluxError << "Target is " << IDcurrentTree << std::endl;
        outputFluxError << "Previous tree was " << IDpreviousTree << std::endl;
        outputFluxError << "Target coordinates are " << agentCurrentLocation[0] << " " << agentCurrentLocation[1] << std::endl;
        outputFluxError << "Previous coordinates were " << agentPreviousLocation[0] << " " << agentPreviousLocation[1] << std::endl;
        outputFluxError << "Distance to travel is " << distanceTravelled << std::endl;
        outputFluxError << "Foraged based on memory " << forageBasedOnMemory << std::endl;
        Rcpp::stop("Problem, negative or null distance");
        outputFluxError.close();
        // Rcout << "Previous target ID: " << IDpreviousTree << " Current target ID: " << IDcurrentTree << std::endl;
      }

      totDistanceTravelled += distanceTravelled;
      // Rcout << "Distance travelled " << distanceTravelled << std::endl;
      // Update time
      timerPrevious = timer;
      timer += distanceTravelled/speed;
    }

    //*****************
    // 2) UPDATE ENVIRONMENT: food depletion/growth
    //*****************

    //a) Update food due to time increase
    //foodConsumedPrevious_v = foodConsumed_v;
    foodQuantityAtTree_v = foodTree(
      timer,
      fruitingTimes_m(_,0),
      fruitingTimes_m(_,1),
      cycleLength,
      foodConsumed_v,
      maximumFoodToYield_v,
      linear
    );

    //b) Update last visit, last visit with feeding and consumed food based on visited pattern (the first two) and time in the year
    NumericMatrix matrixUpdateVisitTimeAndConsumption = updateVisitTimesAndConsumptionCycle(//Triangular distribution of food quantity
      timer,
      fruitingTimes_m(_,0),
      fruitingTimes_m(_,1),
      cycleLength,
      foodConsumed_v,
      visitedTrees_v,
      lastTimeVisitTrees_v,
      lastTimeVisitTreesWithFruit_v,
      foodQuantityAtTree_v
    );

    lastTimeVisitTrees_v = matrixUpdateVisitTimeAndConsumption(_,0);
    // if(timer >= 7206){
    //   Rcout << "Number visited trees at bout " << sum(visitedTrees_v) << std::endl;
    //   Rcout << "Number visited trees in memory " << sum(ifelse(timer - lastTimeVisitTrees_v < noReturnTime, 1, 0)) << std::endl;
    // }
    lastTimeVisitTreesWithFruit_v = matrixUpdateVisitTimeAndConsumption(_,1);
    foodConsumed_v = matrixUpdateVisitTimeAndConsumption(_,2);

    //c) Save eaten food after increase
    //Rcout << "Quantity food " << sum(foodQuantityAtTree_v) << std::endl;
    NumericVector foodEatenDuringBout_v = foodQuantityAtTree_v * as<NumericVector>(visitedTrees_v);
    //Rcout << "Food eaten " << sum(foodEatenDuringBout_v) << std::endl;
    totFoodEaten += sum(foodEatenDuringBout_v);
    //Rcout << "Distance and food tot " << totDistanceTravelled << " " << totFoodEaten << std::endl;

    //d) Update food at tree after depletion
    //Rcout << "Quantity food visited trees prior to update: " << sum(ifelse(visitedTrees_v == 1, foodQuantityAtTree_v, 0)) << std::endl;
    foodQuantityAtTree_v = foodTree(
      timer,
      fruitingTimes_m(_,0),
      fruitingTimes_m(_,1),
      cycleLength,
      foodConsumed_v,
      maximumFoodToYield_v,
      linear
    );

    // if(timer >= 7206){
    //   for(int tree = 0; tree < numberTrees; tree++){
    //     if(visitedTrees_v[tree] == 1){
    //       Rcout << "Was visited " << visitedTrees_v[tree] << " Food presence " << foodQuantityAtTree_v[tree] << std::endl;
    //     }
    //   }
    // }
    //Rcout << "Quantity food visited trees after update: " << sum(ifelse(visitedTrees_v == 1, foodQuantityAtTree_v, 0)) << std::endl;

    if((sum(ifelse(visitedTrees_v == 1, foodQuantityAtTree_v, 0)) > 0.01)){
      IntegerVector vectorTreeWithProblem = rep(0, visitedTrees_v.size());
      for(int vec = 0; vec < vectorTreeWithProblem.size(); vec++){
        if(visitedTrees_v[vec] == 1 && foodQuantityAtTree_v[vec] > 0.01){
          vectorTreeWithProblem[vec] = 1;
        }
      }

      int IDTreeProblem = which_max(vectorTreeWithProblem);
      Rcout << "Problem " << IDTreeProblem << " " << fruitingTimes_m(IDTreeProblem,0) << " " << fruitingTimes_m(IDTreeProblem,1) << " " << foodQuantityAtTree_v[IDTreeProblem] << std::endl;
      Rcout << " " << visitedTrees_v[IDTreeProblem] << " " << foodConsumed_v[IDTreeProblem] << std::endl;
      Rcout << IDpreviousTree << " " << IDcurrentTree << std::endl;
      //I have checked. Calculations are good; it is a very low number when it occurs -> certainly a problem when substraction occurs, residuals > 0.01 kept
    }

   // if(timer >= 7206){
   //  Rcout << "Everything updated but memory" << std::endl;
   //  for(int tree = 0; tree < numberTrees; tree++){
   //    if(visitedTrees_v[tree] == 1){
   //      Rcout << "Was visited " << visitedTrees_v[tree] << " Last time visit " << lastTimeVisitTrees_v[tree] - timer << std::endl;
   //    }
   //  }
   // }

   //*****************
   // 3) UPDATE AGENT KNOWLEDGE
   //*****************

   if(learning){
     //Update memory: visited fruiting trees can be memorised to the detriment of trees visited the longest time ago
     for(int updateMemoryCount = 0; updateMemoryCount < numberTrees; updateMemoryCount++){
      if(spatialKnowledgeRate > 0.01){
       if(lastTimeVisitTreesWithFruit_v[updateMemoryCount] == timer & spatiotemporalKnowledge_m(updateMemoryCount, 0) == 0){//if dead tree not known
         //lastTimeVisitTrees_v[updateMemoryCount]
         //Change the time of last visit of unknown trees not to take them into account for memory lapse
         NumericVector lastTimeVisitTrees_v_memorisedTrees = ifelse(spatiotemporalKnowledge_m(_, 0) == 0, timer + 1, lastTimeVisitTrees_v);
         int treeOldestInMemory = which_minRcpp(lastTimeVisitTrees_v_memorisedTrees - timer)[0];
         // Rcout << treeOldestInMemory << std::endl;
         // Rcout << lastTimeVisitTrees_v[treeOldestInMemory] << std::endl;
         // Rcout << treeLocReal_m(updateMemoryCount, 0) << " " << treeLocReal_m(updateMemoryCount, 1) << std::endl;

         if(lastTimeVisitTreesWithFruit_v[treeOldestInMemory] != timer){
           //lastTimeVisitTreesWithFruit_v[updateMemoryCount]
           //Memorise seedling
           spatiotemporalKnowledge_m(updateMemoryCount, 0) = 1;
           spatiotemporalKnowledge_m(updateMemoryCount, 1) = spatiotemporalKnowledge_m(treeOldestInMemory, 1);//In case you have a tree memorised spatially but not temporally, temporality is only memorised if that was the case of the tree visited the longest time ago. For the moment useless since spatial and temporal memory are equivalent.

           //Erase tree visited the longest time ago
           spatiotemporalKnowledge_m(treeOldestInMemory, 0) = 0;
           spatiotemporalKnowledge_m(treeOldestInMemory, 1) = 0;
         }else{
           Rcout << "Problem in attributing new memorised tree" << std::endl;
         }
       }
      } else{
        //Do nothing
      }
     }
   }

  //*****************
  // 4) DISPERSAL
  //*****************

  if(timer < cycleLimitNumber * cycleLength){//Dispersal only occur in the "true" simulation but not the last run, meant to analyse routine
     if(enterTorpor){

     }else{
       //Re-update if already dispersed binary vector
       disseminated_v = ifelse(lastTimeDispersal_v < timer - timeDelayForDispersal, 0, disseminated_v);

       //Update distribution after seed dispersal
       if(timer > 0.01){

         //Modify if Dispersal realised
         int counterDispersedTrees(0);
         for(int DispersalLoopCount = 0; DispersalLoopCount < numberTrees; DispersalLoopCount++){
           double probaDispersalRealised = runif(1, 0, 1)[0];
           if((lastTimeVisitTreesWithFruit_v[DispersalLoopCount] >= timer - timeDelayForDispersal) & // Tree visited with fruit
              (probaDispersalRealised <= DispersalProbability * (timer - timerPrevious)/timeDelayForDispersal) & //Proba Dispersal realised
              (disseminated_v[DispersalLoopCount] == 0) //Was not already disseminated
           ){
             if(samplingMapCounter==0){
               Rcout << "ERROR: Dispersal occurs at start" << std::endl;
             }
             // if(timer >= 7206){
             //  Rcout << "ID " << DispersalLoopCount << " dispersed" << std::endl;
             // }
             if(distanceTravelled > sqrt(spaceTreeSq)){//No dispersal if movement is insufficient for a free space to be available
             //~~~~~~~~
             // IF RANGE DISPERSAL IS NULL: Dispersal on the path line
             //~~~~~~~~
             // In case dispersal occur only on the path line, find the "available" space using a systemic approach to speed up

             double angleToTarget = atan2((agentCurrentLocation[1] - agentPreviousLocation[1]),(agentCurrentLocation[0] - agentPreviousLocation[0]));
             NumericVector stepTravel_v = rcpp_seq(0, distanceTravelled, 1);//A unitary step of travel to create coordinates along the way.
             // Rcout << distanceTravelled << std::endl;
             // Rcout << "Seq travel " << min(stepTravel_v) << " " << max(stepTravel_v) << std::endl;
             // Rcout << "Angle " << angleToTarget << std::endl;
             // Rcout << "Loc start " << agentPreviousLocation[0] << " " << agentPreviousLocation[1]  << std::endl;
             // Rcout << "Loc end " << agentCurrentLocation[0] << " " <<  agentCurrentLocation[1] << std::endl;

             //Sequence of all coordinates
             NumericVector possibleCoordinatesX_v = agentPreviousLocation[0] + stepTravel_v * cos(angleToTarget);//rcpp_seq(agentPreviousLocation[0], agentCurrentLocation[0], cos(angleToTarget));
             NumericVector possibleCoordinatesY_v = agentPreviousLocation[1] + stepTravel_v * sin(angleToTarget);//rcpp_seq(agentPreviousLocation[1], agentCurrentLocation[1], sin(angleToTarget));

             // Rcout << "X coord " << min(possibleCoordinatesX_v) << " " << max(possibleCoordinatesY_v) << std::endl;

             //The coordinates for the seedling: will only be effective once (and if) free space is found (in that case coordinates are changed)
             bool hasChangedInitialCoordinates(false);
             NumericVector newCoordinatesAfterDispersal_v = rep(0.0, 2);//will only be used if coordinates change

             //Sample only among those that are at sufficient distance of other trees
             NumericVector distanceClosestTree_v = rep(0.0, possibleCoordinatesX_v.size());
             IntegerVector isAboveDistance = rep(0, possibleCoordinatesX_v.size());
             //Rcout << "Worked up to here" << std::endl;
             IntegerVector whichToChoose_v;

             for(int loc = 0; loc < possibleCoordinatesX_v.size(); loc++){
               distanceClosestTree_v[loc] = min(
                 pow(treeLocReal_m(_,0) - possibleCoordinatesX_v[loc],2) +
                   pow(treeLocReal_m(_,1) - possibleCoordinatesY_v[loc],2)
               );
               if(distanceClosestTree_v[loc] > spaceTreeSq){
                 isAboveDistance[loc] = 1;
                 whichToChoose_v.push_back(loc);
               }
             }

             // Rcout << sum(isAboveDistance) << std::endl;
             //Sample one set of coordinates for those that are not too close
             if(sum(isAboveDistance) >= 1){
               int whichLocToChoose = sample(whichToChoose_v, 1)[0];
               newCoordinatesAfterDispersal_v[0] = possibleCoordinatesX_v[whichLocToChoose];
               newCoordinatesAfterDispersal_v[1] = possibleCoordinatesY_v[whichLocToChoose];
               hasChangedInitialCoordinates = true;
               // if(timer >= 7206){
               //   Rcout << "Found dispersal location" << std::endl;
               // }
             }else{
               //No dispersal
               // if(timer >= 7206){
               //  Rcout << "FAILED DISPERSAL" << std::endl;
               // }
             }

               //~~~~~~~~
               // IF RANGE DISPERSAL IS NOT NULL
               //~~~~~~~~
               // //Get the new coordinates:only if space is empty, otherwise sprout does not grow and dispersal not considered
               // int trials(0);
               // NumericVector newCoordinatesAfterDispersal_v = rep(0.0, 2);
               // bool isTreeNearby(true);
               // bool hasChangedInitialCoordinates(false);
               //
               // while(isTreeNearby && trials < 999){//Control if range dispersal is not 0
               //   //Rcout << "Trying to find empty space " << trials << std::endl;
               //   if(trials ==1){
               //     Rcout << "Difficulty to find empty space" << std::endl;
               //   }
               //  if(trials == 998){
               //    Rcout << "Empty space not found" << std::endl;
               //  }
               //   trials += 1;
               //   newCoordinatesAfterDispersal_v = newCoordinatesAfterDispersal(
               //     agentCurrentLocation,
               //     agentPreviousLocation,
               //     0,
               //     mapSize
               //   );
               //
               //   NumericVector distanceClosestTree_v =
               //     (pow(treeLocReal_m(_, 0) - newCoordinatesAfterDispersal_v[0], 2) +
               //     pow(treeLocReal_m(_, 1) - newCoordinatesAfterDispersal_v[1], 2));
               //
               //   int whichMinDistanceTree = which_minRcpp(distanceClosestTree_v)[0];
               //   // Rcout << whichMinDistanceTree << std::endl;
               //   // if(trials == 1){
               //   //   Rcout << "Coord x " << newCoordinatesAfterDispersal_v[0] << std::endl;
               //   //   Rcout << "Coord y " << newCoordinatesAfterDispersal_v[1] << std::endl;
               //   //   Rcout << "Dist closest " << distanceClosestTree_v[whichMinDistanceTree] << std::endl;
               //   // }
               //   if(distanceClosestTree_v[whichMinDistanceTree] > (spaceTreeSq * spaceTreeSq)){
               //     isTreeNearby = false;
               //     hasChangedInitialCoordinates = true;
               //     Rcout << "Distance is " << distanceClosestTree_v[whichMinDistanceTree] << std::endl;
               //     Rcout << "Should be superior to " << spaceTreeSq * spaceTreeSq << std::endl;
               //   }
               // }

               if(hasChangedInitialCoordinates){

                 //If the trees was previously visited, with fruits, and seed "effectively" digested/transported and germinated (i.e. probability of Dispersal realised), then change coordinates:

                 //^^^^^^^^^^
                 // IN CASE NEW TREE REMPLACES FATHER/MOTHER TREE
                 //^^^^^^^^^^

                 // treeLocReal_m(DispersalLoopCount, 0) = newCoordinatesAfterDispersal_m(DispersalLoopCount, 0);
                 // treeLocReal_m(DispersalLoopCount, 1) = newCoordinatesAfterDispersal_m(DispersalLoopCount, 1);
                 // //Change if was disseminated to avoid redissemination
                 // disseminated_v[DispersalLoopCount] = 1;
                 // lastTimeDispersal_v[DispersalLoopCount] = timer;
                 // lastTimeVisitTrees_v[DispersalLoopCount] = timer;

                 //^^^^^^^^^^
                 // IN CASE NEW TREE REMPLACES RANDOMLY ANOTHER TREE
                 //^^^^^^^^^^

                 //Randomly chose tree that dies and is replaced by new one
                 int treeToDelete = floor(runif(1, 0, numberTrees)[0]);
                 while(treeToDelete == IDcurrentTree | treeToDelete == IDpreviousTree){
                   treeToDelete = floor(runif(1, 0, numberTrees)[0]);
                 }
                 //Rcout << "I changed tree number: " << treeToDelete << std::endl;
                 treeLocReal_m(treeToDelete, 0) = newCoordinatesAfterDispersal_v[0];
                 treeLocReal_m(treeToDelete, 1) = newCoordinatesAfterDispersal_v[1];

                 //Resample start date according to father/mother one
                 double newStartDate = rnorm(1, fruitingTimes_m(DispersalLoopCount, 0), 1)[0];
                 if(newStartDate < 0){
                   newStartDate = fmod(newStartDate + cycleLength*10, cycleLength);
                 }
                 Rcout << "Mother tree fruiting date " << fruitingTimes_m(DispersalLoopCount, 0) << std::endl;
                 Rcout << "New tree fruiting date " << newStartDate << std::endl;
                 fruitingTimes_m(treeToDelete, 0) = newStartDate;
                 fruitingTimes_m(treeToDelete, 1) = newStartDate + fruitingLength;

                 //Change if was disseminated to avoid redissemination
                 disseminated_v[DispersalLoopCount] = 1;
                 lastTimeDispersal_v[DispersalLoopCount] = timer;
                 lastTimeVisitTrees_v[DispersalLoopCount] = timer;

                 disseminated_v[treeToDelete] = 1;
                 lastTimeDispersal_v[treeToDelete] = timer;
                 lastTimeVisitTrees_v[treeToDelete] = timer;

                 //Fake food consumption to avoid visit based on memory if new date is posterior to current time
                 foodConsumed_v[treeToDelete] = 1;

                 if(learning){
                   //Update memory
                   //If the now "dead" tree was memorised, so is the seedling. Otherwise, the tree visited the longest time ago will be. The subsequent lines deal with this case.
                   if(spatiotemporalKnowledge_m(treeToDelete, 0) == 0){//if dead tree not known

                     //Change the time of last visit of unknown trees not to take them into account for memory lapse
                     NumericVector lastTimeVisitTrees_v_memorisedTrees = ifelse(spatiotemporalKnowledge_m(_, 0) == 0, timer + 1, lastTimeVisitTrees_v);
                     int treeOldestInMemory = which_minRcpp(lastTimeVisitTrees_v_memorisedTrees - timer)[0];

                     //Memorise seedling
                     spatiotemporalKnowledge_m(treeToDelete, 0) = 1;
                     spatiotemporalKnowledge_m(treeToDelete, 1) = spatiotemporalKnowledge_m(treeOldestInMemory, 1);//In case you have a tree memorised spatially but not temporally, temporality is only memorised if that was the case of the tree visited the longest time ago. For the moment useless since spatial and temporal memory are equivalent.

                     //Erase tree visited the longest time ago
                     spatiotemporalKnowledge_m(treeOldestInMemory, 0) = 0;
                     spatiotemporalKnowledge_m(treeOldestInMemory, 1) = 0;
                   }
                 }
                 counterDispersedTrees += 1;
               }
             }
           }
           //Rcout << counterDispersedTrees << " trees have dispersed!" << std::endl;
         }
       }else{

       }
       // if(timer >= 7206){
       //  Rcout << "Everything updated" << std::endl;
       // }
       // if(timer >= 7206){
       //   //Check whether a) number of known trees is kept constant b) spatial = temporal knowledge. If so, it means that replacement above, and due to seed dispersal, is correctly done.
       //   Rcout << "Memorised trees " << sum(spatiotemporalKnowledge_m(_, 0)) << " " << sum(spatiotemporalKnowledge_m(_, 1)) << std::endl;
       //
       //   //Chek if IDs known trees change: should indicate new memorisation
       //   int sumIDKnown = 0;
       //   for(int counter = 0; counter < numberTrees; counter++){
       //     if(spatiotemporalKnowledge_m(counter, 0)==1){
       //       sumIDKnown += 1;
       //     }
       //  }
       //   Rcout << "Trees known sum ID " << sumIDKnown << std::endl;
       // }
     }
   }

  ////------
  // Saving efficiency results: semi continuous, each cycle length / 5 tu, up to end of normal run (i.e. when there is dispersal)
  ////------

  if(outputFluxEfficiencyContinuous){
    // NumericVector nearestNeighbourMeanSd = nearestNeighbourDistance(
    //   treeLocReal_m,
    //   mapSize + 1
    // );
    //
    // NumericVector lloydIndices = lloydIndex(
    //   treeLocReal_m,
    //   mapSize,
    //   quadratSize
    // );
    //
    // Rcout << "Calculate spatial metrics" << std::endl;

    //Rcout << timer << " " << fmod(timer, cycleLength/5) << " " << fmod(timerPrevious, cycleLength/5) << std::endl;
    if(fmod(timer, cycleLength/5) < fmod(timerPrevious, cycleLength/5)){
      //Rcout << "I should save" << std::endl;
      outputFluxEfficiencyContinuous <<
        timer << " " <<
          totFoodEaten/totDistanceTravelled << std::endl;//" " <<
            // nearestNeighbourMeanSd[0] << " " <<
            //   nearestNeighbourMeanSd[1] << " " <<
            //     lloydIndices[0] << " " <<
            //       lloydIndices[1] << std::endl;
    }
    // if(timer >= 7206){
    //   Rcout << "Continuous output updated" << std::endl;
    // }
  }

  //----------------
  //Save results map
  //----------------

  if(saveTreeMap){
  if(
    (samplingMapCounter == 0) |
    ((floor(timer) >= samplingMapTime_v[samplingMapCounter]) & (floor(timerPrevious) <= samplingMapTime_v[samplingMapCounter]))
  ){
    // Rcout << "Saving map " << samplingMapTime_v[samplingMapCounter] << std::endl;
    // Rcout << timer << std::endl;
    //Update counter to avoid resampling map multiple times
    samplingMapCounter += 1;

    ////------
    // Saving tree map
    ////------

    std::string pathOfTheFileOutputInitMap = pathOfTheFileOutputInit;
    pathOfTheFileOutputInitMap.append("_Map_");
    pathOfTheFileOutputInitMap.append(std::to_string(samplingMapTime_v[samplingMapCounter - 1]));
    pathOfTheFileOutputInitMap.append(".txt");

    std::ofstream outputFluxMap;
    //Open connection to output file
    outputFluxMap.open(pathOfTheFileOutputInitMap.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
    if(outputFluxMap)  // Write the file only if correctly opened
    {
      //Column names
      outputFluxMap <<
        "x" << " " <<
          "y" << " " <<
            "startFruit" << " " <<
              "endFruit" << std::endl;

      //Values
      for(int loopingForMapSave = 0; loopingForMapSave < numberTrees; loopingForMapSave++){
        //Rcout << loopingForMapSave << std::endl;
        outputFluxMap <<
          treeLocReal_m(loopingForMapSave, 0) << " " <<
            treeLocReal_m(loopingForMapSave, 1) << " " <<
              fruitingTimes_m(loopingForMapSave,0) << " " <<
                fruitingTimes_m(loopingForMapSave,1) << std::endl;
      }
    }
    else //If impossible to write in the file, say it
    {
      Rcout << "ERROR: Impossible to open the output file about mapping" << std::endl;
    }

    //Close connection to output file
    outputFluxMap.close();
    }
  }

  // Rcout << "End bout" << std::endl;

  ////------
  // Saving efficiency results: end only
  ////------

  if(timer >= cycleLimitNumber * cycleLength && hasSavedFinalResults == 0){
  hasSavedFinalResults = 1;
  //Close connection to output file
  outputFluxEfficiencyContinuous.close();

  std::string pathOfTheFileOutputInitEfficiency = pathOfTheFileOutputInit;
  pathOfTheFileOutputInitEfficiency.append("_Efficiency.txt");

  std::ofstream outputFluxEfficiency;
  //Open connection to output file
  outputFluxEfficiency.open(pathOfTheFileOutputInitEfficiency.c_str(), std::ios::out | std::ios::trunc);//c_str() convert it to a string category understandable for the ofstream function, first creates file or blanks existing one
  if(outputFluxEfficiency)  // Write the file only if correctly opened
  {
    //Column names
    outputFluxEfficiency <<
      "Run_simulation" << " " <<
        "Length_run" << " " <<
          "Cycle_length" << " " <<
            "Fruiting_length" << " " <<
              "Map_size" << " " <<
                "Ntrees" << " " <<
                  "Distribution_homogeneous_at_start" << " " <<
                    "Nclusters_trees" << " " <<
                      "Spread_clusters_trees" << " " <<
                        "Maxyieldfood" << " " <<
                          "Speed" << " " <<
                            "Torpor_duration" << " " <<
                              "No_return_time" << " " <<
                                "Value_if_unknown_temporality" << " " <<
                                  "What_rule_to_move" << " " <<
                                    "intensityCompetitionForSpace" << " " <<
                                      "moveOnlyToFruitingTrees" << " " <<
                                        "moveOnlyToTarget" << " " <<
                                          "linear" << " " <<
                                            "learning" << " " <<
                                              "Dispersal_probability" << " " <<
                                                "Time_delay_Dispersal" << " " <<
                                                  "Perception_range" << " " <<
                                                    "Temporal_knowledge_rate" << " " <<
                                                      "Spatial_knoledge_rate" << " " <<
                                                        "Number_events_with_no_food" << " " <<
                                                          "Tot_food_eaten" << " " <<
                                                            "Tot_distance_travelled" << std::endl;

    //Values
    outputFluxEfficiency <<
      repetitionNumber << " " <<
        timer << " " <<
          cycleLength << " " <<
            fruitingLength << " " <<
              mapSize << " " <<
                numberTrees << " " <<
                  homogeneousDistribution << " " <<
                    treeClusterNumber << " " <<
                      treeClusterSpread << " " <<
                        max(maximumFoodToYield_v) << " " <<
                          speed << " " <<
                            torporTime << " " <<
                              noReturnTime << " " <<
                                whatValueUnknownTemporal << " " <<
                                  whatRule << " " <<
                                    intensityCompetitionForSpace << " " <<
                                      moveOnlyToFruitingTrees << " " <<
                                        moveOnlyToTarget << " " <<
                                          linear << " " <<
                                            learning << " " <<
                                              DispersalProbability << " " <<
                                                timeDelayForDispersal << " " <<
                                                  perceptualRange << " " <<
                                                    temporalKnowledgeRate << " " <<
                                                      spatialKnowledgeRate << " " <<
                                                        eventsWithNoFood << " " <<
                                                          totFoodEaten << " " <<
                                                            totDistanceTravelled << std::endl;

  }
  else //If impossible to write in the file, say it
  {
    Rcout << "ERROR: Impossible to open the output file about foraging efficiency" << std::endl;
  }

  //Close connection to output file
  outputFluxEfficiency.close();
  }
  }
  //Close connection to output file routine
  outputRoutine.close();
}