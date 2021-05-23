#include "Tangle.h"
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <numeric>
#include <fstream>
#include <cmath>
#include <functional>

int TxActor::actorCount;
long int Tx::tx_totalCount;
int Tangle::TangleGiveTipsCount;

/*
    Tx DEFINITIONS
*/

//return true if Tx has approved at least one other transactions
bool Tx::hasApprovees()
{
    return m_TxApproved.size() > 0;
}

Tx::Tx() : m_walkBacktracks(0), TxNumber(tx_totalCount)
{
    tx_totalCount++;
}


//Tx def END

/*
    TANGLE DEFINITIONS
*/

Tangle::Tangle() try : m_genesisBlock( new Tx )
{
     m_tips[m_genesisBlock->TxNumber] =  m_genesisBlock;
     unsigned tipSelectSeed;
     tipSelectSeed = std::chrono::system_clock::now().time_since_epoch().count(); //needs to use seeds from omnetpp
     tipSelectGen.seed( tipSelectSeed );
     m_genesisBlock->isGenesisBlock = true;
}

catch ( std::bad_alloc e)
{
     std::cerr << e.what() << std::endl;
     std::cerr << "Failed to create genesis transaction" << std::endl;
}

//method to return a copy of all the current unconfirmed transactions
std::map<int, t_ptrTx> Tangle::giveTips()
{
     return m_tips;
}

//checks the newly approved tips against the current tip vecotr held by tangle
//removes any that are still unconfirmed
void Tangle::ReconcileTips( const t_txApproved& removeTips )
{

    for(auto& tipSelected : removeTips)
    {

        auto it = m_tips.find( tipSelected->TxNumber );

        if( it != m_tips.end() )
        {
            m_tips.erase( it );
        }
    }

}

//adds a pointer to a newly added but as yet unconfirmed Tx to the tip list
void Tangle::addTip( t_ptrTx newTip )
{
     m_tips[newTip->TxNumber] = newTip;
     allTx.push_back( newTip );
}

std::mt19937& Tangle::getRandGen()
{
    return tipSelectGen;
}

int Tangle::getTipNumber()
{
    return m_tips.size();
}

const t_ptrTx& Tangle::giveGenBlock() const
{
    return m_genesisBlock;
}

void Tangle::printTangle(std::vector<t_ptrTx> allTx,const t_ptrTx& Gen)
{
    std::fstream file;
    std::string path = "./data/Tracking/TrackerTangle.txt";
    file.open(path,std::ios::app);

    omnetpp::simtime_t sec = omnetpp::simTime();

    file << sec << ";";
    file <<"Genesis"<< ";";

    for(long unsigned int i = 0; i < Gen->m_approvedBy.size(); i++)
    {
        if(i == Gen->m_approvedBy.size() - 1)
        {
            file << Gen->m_approvedBy[i]->TxNumber + 1;
        }

        else
        {
            file << Gen->m_approvedBy[i]->TxNumber + 1 << ",";
        }
    }

    file << std::endl;

    for(long unsigned int i = 0; i < allTx.size(); i++)
    {
       file << sec << ";";
       file << allTx[i]->TxNumber + 1 << ";";

        for(long unsigned int j = 0; j < allTx[i]->m_approvedBy.size(); j++)
        {
            auto temp = allTx[i]->m_approvedBy[j];

            if(j == allTx[i]->m_approvedBy.size() - 1)
            {
                file << temp->TxNumber + 1;
            }

            else
            {
                file << temp->TxNumber + 1 << ",";
            }
        }

        file << std::endl;
    }

    file << std::endl;

    file.close();
}

void Tangle::printTipsLeft(int numberTips)
{
    std::fstream file;
    std::string path = "./data/log/NumberTips.txt";
    file.open(path,std::ios::app);

    file << numberTips << std::endl;

    file.close();
    return;
}

// Tangle def END

/*
    TXACTOR DEFINITIONS
*/

//Constructor

TxActor::TxActor() { ++actorCount; }

//Tips to approve selected completely at random
t_txApproved TxActor::URTipSelection( std::map<int, t_ptrTx> tips )
{

     t_txApproved chosenTips;

     for ( int i = 0; i < APPROVE_VAL; ++i )
     {
         if(tips.size() > 0)
         {
             std::uniform_int_distribution<int> tipDist( 0, tips.size() -1 );
             int iterAdvances = tipDist( getTanglePtr()->getRandGen() );

             assert(iterAdvances < tips.size());

             auto beginIter = tips.begin();
             std::advance( beginIter, iterAdvances );

             chosenTips.push_back( beginIter->second );
             tips.erase( beginIter );

             if( tips.size() == 0 )
             {
                 break;
             }

         }
     }

     chosenTips.erase( std::unique( chosenTips.begin(), chosenTips.end() ), chosenTips.end() ) ;
     return chosenTips;
}

//creates a new transaction, selects tips for it to approve, then adds the new transaction to the tip list
//ready for approval by the proceeding transactions
void TxActor::attach( std::map<int, t_ptrTx>& storedTips, omnetpp::simtime_t attachTime, t_txApproved& chosen )
{
     try
     {

         //create new tx
         m_MyTx.emplace_back( new Tx() );
         m_MyTx.back()->m_issuedBy = this;
         m_MyTx.back()->timeStamp = attachTime;


         //add pointer to new Tx to tips selected, so they know who approved them
         for ( auto& tipSelected : chosen )
         {

             tipSelected->m_approvedBy.push_back( m_MyTx.back() );

             if( !( tipSelected->isApproved ) )
             {
                 //with firstApprovedTime and timeAttached as field - we can compute the age of a transaction
                 tipSelected->firstApprovedTime = attachTime;
                 tipSelected->isApproved = true;
             }

         }

         m_MyTx.back()->m_TxApproved = chosen;

         //remove pointers to tips just approved, from tips vector in tangle
         getTanglePtr()->ReconcileTips( chosen );

         //add newly created Tx to Tangle tips list
         getTanglePtr()->addTip( m_MyTx.back() );

     }
     catch ( std::bad_alloc& e )
     {
         std::cerr << e.what() << std::endl;
         std::cerr << "Not enough memory for TxActor to produce a new Tx" << std::endl;
     }

}

//getter and setter for tangle pointer, makes it easy to interact with the tips
 Tangle* TxActor::getTanglePtr() const
{
    return tanglePtr;
}

void TxActor::setTanglePtr( Tangle* tn )
{
    tanglePtr = tn;
}

const std::vector<t_ptrTx>& TxActor::getMyTx() const
{
    return m_MyTx;
}

//compute weight definitions

//indirect recursion, start point here - calls private recursive function _computeWeight
int TxActor::ComputeWeight( t_ptrTx tx, omnetpp::simtime_t timeStamp )
{

    std::vector<t_ptrTx> visited;
    int weight = _computeWeight( visited, tx, timeStamp );

    //leave txes as we found them
    for( int i = 0; i < visited.size(); ++i )
    {

        for( int j = 0; j < visited.at( i )->m_approvedBy.size(); ++j )
        {
            visited.at( i )->m_approvedBy.at( j )->isVisited = false;
        }

    }

    tx->isVisited = false;
    return weight + 1;

}

//traverses tangle returns weight of each transaction, stopping cases: previously visited transaction, transaction with
//a timestamp after TxActor started computing, and on reaching a tip
int TxActor::_computeWeight( std::vector<t_ptrTx>& visited, t_ptrTx& current, omnetpp::simtime_t timeStamp )
{

    //could TxActor "see" the current Tx their walk is located on?
    if( timeStamp < current->timeStamp )
    {
        visited.push_back( current );
        current->isVisited = true;
        return 0;
    }

    //if we reach a tip
    if( current->m_approvedBy.size() == 0 )
    {
        visited.push_back( current );
        current->isVisited = true;

        return 0;
    }

    visited.push_back( current );

    current->isVisited = true;
    int weight = 0;

    for( int i = 0; i < current->m_approvedBy.size(); ++i )
    {

        if( !current->m_approvedBy.at( i )->isVisited )
        {
                //check if next tx has been visited before
                weight += 1 + _computeWeight( visited, current->m_approvedBy.at( i ), timeStamp );
        }

    }

    return weight;

}

t_ptrTx TxActor::getWalkStart( std::map<int, t_ptrTx>& tips, int backTrackDist)
{
    std::uniform_int_distribution<int> tipDist( 0, tips.size() -1 );
    int iterAdvances = tipDist( getTanglePtr()->getRandGen() );

    assert( iterAdvances < tips.size() );

    auto beginIter = tips.begin();

    if(tips.size() > 1)
    {
        std::advance( beginIter, iterAdvances );
    }

    t_ptrTx current = beginIter->second;

    int count = backTrackDist;

    int approvesIndex;

    //start backtrack
    //go until genesis block or reach backtrack distance
    while( !current->isGenesisBlock && count > 0 )
    {

        std::uniform_int_distribution<int> choice( 0, current->m_TxApproved.size() -1 );
        approvesIndex = choice( getTanglePtr()->getRandGen() ) ;

        assert( approvesIndex < current->m_TxApproved.size() );

        current = current->m_TxApproved.at( approvesIndex );


        --count;

    }
    
    return current;
}

t_ptrTx TxActor::WalkTipSelection( t_ptrTx start, double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp )
{

    // Used to determine the next Tx to walk to
    int walkCounts = 0;

    t_ptrTx current = start;

    //keep going until we reach a "tip" in relation to the view of the tangle that TxActor has
    while( !isRelativeTip( current, tips ) )
    {
        ++walkCounts;

        //copy of each transactions approvers
        std::vector<t_ptrTx> currentView = current->m_approvedBy;

        //filter current View
        filterView( currentView, timeStamp );

        if( currentView.size() == 0 )
        {
            break;
        }

        //if only one approver available dont compute the weight
        if( currentView.size() == 1 )
        {
            current = currentView.front();

        }
        else
        {
            //if more than one find the max weight
            //get heaviest tx to simplify choosing next
            int maxWeightIndex = findMaxWeightIndex( currentView, timeStamp );
            t_ptrTx heaviestTx = currentView.at( maxWeightIndex );

            currentView.erase( currentView.begin() + maxWeightIndex );

            //if no more after heaviest removal, choose heaviest
            if( currentView.size() == 0 )
            {
                current = heaviestTx;

            }
            else
            {
                // if at least one non heaviest still available pick between them
                std::uniform_real_distribution<double> walkChoice( 0.0, 1.0 );

                if( walkChoice( getTanglePtr()->getRandGen() ) < alphaVal)
                {
                    current = heaviestTx;
                }
                else
                {
                    //otherwise choose a random site
                    if( currentView.size() == 1 )
                    {
                        //if only one left after removing heaviest use the remaining one
                        current = currentView.front();

                    }
                    else
                    {
                        //otherwise pick at random
                        std::uniform_int_distribution<int> siteChoice( 0, currentView.size() - 1 );
                        int choiceIndex = siteChoice( getTanglePtr()->getRandGen() );
                        current = currentView.at( choiceIndex );

                    }
                }
            }
        }
    }

    return current;
}

bool TxActor::isRelativeTip( t_ptrTx& toCheck, std::map<int, t_ptrTx>& tips )
{
    // if tips is sorted when txactor receives it, we can improve performance as this is called alot during tip selection
    auto it = tips.find( toCheck->TxNumber );

    if( it == tips.end() )
    {
        return false;
    }
    else
    {
        return true;
    }

}

void TxActor::filterView( std::vector<t_ptrTx>& view, omnetpp::simtime_t timeStamp )
{

    std::vector<int> removeIndexes;

    for( int i = 0; i < view.size(); ++i )
    {
        if( view.at( i )->timeStamp > timeStamp )
        {
            removeIndexes.push_back( i );
        }
    }

    if( removeIndexes.size() > 0 )
    {
        for( int i = removeIndexes.size() -1; i > -1; --i )
        {
            view.erase( view.begin() + removeIndexes.at( i ) );
        }
    }

}

int TxActor::findMaxWeightIndex(std::vector<t_ptrTx>& view, omnetpp::simtime_t timeStamp )
{
    int maxWeight = 0;
    int maxWeightIndex = 0;

    for( int i = 0; i < view.size(); ++i )
    {
        int weight = ComputeWeight( view.at( i ), timeStamp );

        if( weight > maxWeight )
        {
            maxWeightIndex = i;
        }

    }

    return maxWeightIndex;
}

t_ptrTx TxActor::EasyWalkTipSelection( t_ptrTx start, double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp )
{

    // Used to determine the next Tx to walk to
    int walkCounts = 0;

    t_ptrTx current = start;

    //keep going until we reach a "tip" in relation to the view of the tangle that TxActor has
    while( !isRelativeTip( current, tips ) )
    {

        ++walkCounts;

        //copy of each transactions approvers
        std::vector<t_ptrTx> currentView = current->m_approvedBy;

        //filter current View
        filterView( currentView, timeStamp );

        if( currentView.size() == 0 )
        {
            break;
        }

        //if only one approver available dont compute the weight
        if( currentView.size() == 1 )
        {
            current = currentView.front();

        }
        else
        {

        // choose if calculate the weights of available transaction here
        std::uniform_real_distribution<double> walkChoice( 0.0, 1.0 );

            if( walkChoice( getTanglePtr()->getRandGen() ) < alphaVal)
            {
                //if more than one find the max weight
                //get heaviest tx to simplify choosing next
                int maxWeightIndex = findMaxWeightIndex( currentView, timeStamp );
                current =  currentView.at( maxWeightIndex );


            }
            else
            {

                //otherwise pick at random
                std::uniform_int_distribution<int> siteChoice( 0, currentView.size() - 1 );
                int choiceIndex = siteChoice( getTanglePtr()->getRandGen() );
                current = currentView.at( choiceIndex );

            }
        }
    }

    current->m_walkBacktracks = walkCounts;

    return current;

}

// Allows us to use walk tip selection with multiple walkers
t_txApproved TxActor::NKWalkTipSelection( double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int kMultiplier, int backTrackDist)
{
    // Walkers sent == 3 * k + 4
    int walkers = kMultiplier * APPROVE_VAL + 4;

    std::vector<t_ptrTx> vec_walkerResults;
    vec_walkerResults.reserve(walkers);

    // Let the walkers find the tips tips
    for( int i = 0; i < walkers; ++i )
    {
        vec_walkerResults.push_back( EasyWalkTipSelection( getWalkStart( tips, backTrackDist ) , alphaVal, tips, timeStamp ) );
    }

    // Sort tips by how many steps the walker made - ascending order
    std::sort( vec_walkerResults.begin(), vec_walkerResults.end(), [] ( t_ptrTx left, t_ptrTx right )
        {
            return left->m_walkBacktracks > right->m_walkBacktracks;
        }
    );

    // Dedupe so we dont approve the same tip more than once
    vec_walkerResults.erase( std::unique( vec_walkerResults.begin(), vec_walkerResults.end() ), vec_walkerResults.end() ) ;

    assert( vec_walkerResults.size() > 0 );

    int tipSelectedSize =  APPROVE_VAL > vec_walkerResults.size() ? vec_walkerResults.size() : APPROVE_VAL ;

    t_txApproved retVal(vec_walkerResults.begin(), vec_walkerResults.begin() + tipSelectedSize );

    return retVal;
}

/************************************************ORIGINAL TSA WITH TRACKING**************************************************************/

//backtrack for the WP TSA with the tracking
t_ptrTx TxActor::getWalkStartWPTrack(std::map<int, t_ptrTx>& tips, int backTrackDist, std::string TxActorId)
{
    std::fstream file;
    std::string path = "./data/Tracking/TSATracker" + TxActorId + ".txt";
    file.open(path,std::ios::app);

    std::uniform_int_distribution<int> tipDist( 0, tips.size() -1 );
    int iterAdvances = tipDist( getTanglePtr()->getRandGen() );

    assert( iterAdvances < tips.size() );

    auto beginIter = tips.begin();

    if(tips.size() > 1)
    {
        std::advance( beginIter, iterAdvances );
    }

    t_ptrTx current = beginIter->second;

    if(current->isGenesisBlock)
    {
        file << "   - Selected tip for the back walk: " << "Genesis";
        file.close();
    }

    else
    {
        file << "   - Selected tip for the back walk: " << current->TxNumber + 1;
        file.close();
    }

    int count = backTrackDist;

    int approvesIndex;

    //start backtrack
    //go until genesis block or reach backtrack distance
    while( !current->isGenesisBlock && count > 0 )
    {

        std::uniform_int_distribution<int> choice( 0, current->m_TxApproved.size() -1 );
        approvesIndex = choice( getTanglePtr()->getRandGen() ) ;

        assert( approvesIndex < current->m_TxApproved.size() );

        current = current->m_TxApproved.at( approvesIndex );


        --count;

    }

    return current;
}

//perform the TSA based on the original white paper of Serguei Popov with the tracking
t_ptrTx TxActor::WPWalkTipSelectionTrack(t_ptrTx start, double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int &walk_time, std::string TxActorId)
{
    std::fstream file;
    std::string path = "./data/Tracking/TSATracker" + TxActorId + ".txt";
    file.open(path,std::ios::app);

    // Used to determine the next Tx to walk to
    int walkCounts = 0;
    t_ptrTx current = start;

    if(current->isGenesisBlock)
    {
        file <<"   - Site root: "<< "Genesis" << std::endl;
    }

    else
    {
        file <<"   - Site root: "<< current->TxNumber + 1 << std::endl;
    }

    //keep going until we reach a "tip" in relation to the view of the tangle that TxActor has
    while(!isRelativeTip(current,tips))
    {
        walkCounts++;

        //the weight of current
        int start_weight = ComputeWeight(current,timeStamp);

        if(current->isGenesisBlock)
        {
            file << "       > Current site: " << "Genesis" << ", weight: " << start_weight << ", walk counts: " << walkCounts - 1 << std::endl;
        }

        else
        {
            file << "       > Current site: " << current->TxNumber + 1 << ", weight: " << start_weight << ", walk counts: " << walkCounts - 1 << std::endl;
        }

        file << "       > Adjacent site(s): " << std::endl;

        //copy of each transactions approvers
        std::vector<t_ptrTx> currentView = current->m_approvedBy;

        //store the weight of the sites that approuve current
        std::vector<int> sitesWeight;

        //store the probability to walk in each site that approuve current
        //pair = index,prob
        std::vector<std::pair<int,double>> sitesProb;

        //filter current View
        filterView(currentView,timeStamp);

        if(currentView.size() == 0)
        {
            file << "           -> No adjacent sites available after the call of filterView()" << std::endl;
            file << std::endl;
            break;
        }

        //if only one approver available dont compute the weight
        if(currentView.size() == 1)
        {
            current = currentView.front();
            file << "           -> Only one adjacent site available after the call of filterView():" << std::endl;
            file << "               Site: " << current->TxNumber + 1 << std::endl;
            file << std::endl;
        }

        else
        {
            //choose if calculate the weights of available transaction here
            //random generator
            std::uniform_real_distribution<double> walkChoice(0.0,1.0);

            //the exp sum (see proba formula page 21)
            double sum_exp = 0.0;
            int weight;

            for(int j = 0; j < static_cast<int>(currentView.size()); j++)
            {
                weight = ComputeWeight(currentView[j],timeStamp);
                sum_exp = sum_exp + double(exp(double(-alphaVal*(start_weight - weight))));
                sitesWeight.push_back(weight);
            }

            //a probability (see formula page 21)
            double prob;

            for(int j = 0; j < static_cast<int>(sitesWeight.size()); j++)
            {
               prob = double(exp(double(-alphaVal*(start_weight - sitesWeight[j]))));
               prob = prob/sum_exp;
               sitesProb.push_back(std::make_pair(j,prob));

               file << "           -> Site: " << currentView[j]->TxNumber + 1 <<", weight: " << sitesWeight[j] << ", probability: " << prob << std::endl;
            }

            std::sort(sitesProb.begin(), sitesProb.end(),[](const std::pair<int,double> &a, const std::pair<int,double> &b){
            return a.second < b.second;}); //we sort sitesProb
            int nextCurrentIndex = 0; //index of the selected sites in CurrentView
            double probWalkChoice =  walkChoice(getTanglePtr()->getRandGen()); //the prob of the choosen site for the next walk

            for(int j = 0; j < static_cast<int>(sitesProb.size()) - 1; j++)
            {
                if(j == 0)
                {
                    if(probWalkChoice < sitesProb[j].second)
                    {
                        nextCurrentIndex = sitesProb[j].first;
                        break;
                    }
                }

                else
                {
                    if(probWalkChoice < sitesProb[j-1].second + sitesProb[j].second)
                    {
                        nextCurrentIndex = sitesProb[j].first;
                        break;
                    }
                }
            }

            current = currentView[nextCurrentIndex];
            file << "       > Uniform probability: " << probWalkChoice <<", selected site for the next walk => " << current->TxNumber + 1 << std::endl;
            file << std::endl;
        }
    }



    walk_time = walkCounts;
    current->m_walkBacktracks = walkCounts;

    if(current->isGenesisBlock)
    {
        file << "   - Selected tip for potential approval: " << "Genesis" <<", total walk count:" << walk_time << std::endl;
        file.close();
    }

    else
    {
        file << "   - Selected tip for potential approval: " << current->TxNumber + 1 <<", total walk count: " << walk_time <<std::endl;
        file.close();
    }

    return current;
}

//initiate the TSA with the tracking
t_txApproved TxActor::WhitePaperTSATrack(double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int W, int N, std::string TxActorId)
{
    std::fstream file;
    std::string path = "./data/Tracking/TSATracker" + TxActorId  + ".txt";
    file.open(path,std::ios::app);

    omnetpp::simtime_t sec = omnetpp::simTime();

    file << "\n";
    file <<"* Simulation time: "<< sec <<" sec" << std::endl;
    file << "\n";
    file <<"* Selected sites for the Walk TSA:" << std::endl;

    //random generator with uniform distribution between W and 2*W for sites distance from tips
    std::uniform_int_distribution<int> backTrackDist(W,2*W);

    std::vector<t_ptrTx> start_sites; //vect of start sites for the TSA
    t_ptrTx temp_site;
    int backTD;

    //we pick randomly N sites for the TSA start point (may be the same)
    for(int i = 0; i < N; i++)
    {
        backTD = backTrackDist(getTanglePtr()->getRandGen());
        temp_site = getWalkStartWPTrack(tips,backTrackDist(getTanglePtr()->getRandGen()),TxActorId);
        file <<", back track distance: " << backTD;

        if(temp_site->isGenesisBlock)
        {
            file <<", selected site: " << "Genesis" << std::endl;
            start_sites.push_back(temp_site);
        }

        else
        {
            file <<", selected site: " << temp_site->TxNumber + 1 << std::endl;
            start_sites.push_back(temp_site);
        }

    }

    file << "\n";

    t_txApproved selected_tips; //tips selected by the TSA (may be the same)
    std::vector<std::pair<int,int>> walk_total; //Execution time (i.e walk time) of the TSA for each start sites
    int walk_time; //Execution time for one start site

    file << "* TSA Path:" << std::endl;

    //initiate for each start site a TSA and we save the execution time in time_exe
    for(int i = 0; i < N; i++)
    {
        if(i != N-1)
        {
            selected_tips.push_back(WPWalkTipSelectionTrack(start_sites[i],alphaVal,tips,timeStamp,walk_time,TxActorId));
            file << std::endl;
            walk_total.push_back(std::make_pair(i,walk_time));
        }

        else
        {
            selected_tips.push_back(WPWalkTipSelectionTrack(start_sites[i],alphaVal,tips,timeStamp,walk_time,TxActorId));
            walk_total.push_back(std::make_pair(i,walk_time));
        }
    }

    std::sort(walk_total.begin(), walk_total.end(),[](const std::pair<int,int> &a, const std::pair<int,int> &b){
    return a.second < b.second;});

    t_txApproved tipstoApprove; //the two first tips selected by the TSA based on time execution
    //index in the time_exe vector of the first tip to be selected based on time execution
    int Index_first_tips = walk_total[0].first;
    tipstoApprove.push_back(selected_tips[Index_first_tips]); //we push it

   for(long unsigned int j = 1; j < walk_total.size(); j++)
   {
       Index_first_tips = walk_total[j].first;

       if(selected_tips[Index_first_tips]->TxNumber != tipstoApprove[0]->TxNumber)
       {
           tipstoApprove.push_back(selected_tips[Index_first_tips]);
           break;
       }
   }


    file << "\n";
    file << "* Selected tip(s) for approval:" << std::endl;

    for(long unsigned int k = 0; k < tipstoApprove.size(); k++)
    {
        if(tipstoApprove[k]->isGenesisBlock)
        {
          file <<"    - Tip: " << "Genesis" << std::endl;
        }

        else
        {
          file <<"    - Tip: " << tipstoApprove[k]->TxNumber + 1 << std::endl;
        }
    }

    file << "\n";
    file << "------------------------------------------------------------------------------------";
    file << "\n";

    file.close();

    return tipstoApprove;
}

/************************************************ORIGINAL TSA WITHOUT TRACKING************************************************************/

//perform the TSA based on the original white paper of Serguei Popov without the tracking
t_ptrTx TxActor::WPWalkTipSelection(t_ptrTx start, double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int &walk_time)
{
    // Used to determine the next Tx to walk to
    int walkCounts = 0;
    t_ptrTx current = start;

    //keep going until we reach a "tip" in relation to the view of the tangle that TxActor has
    while(!isRelativeTip(current,tips))
    {
        walkCounts++;

        //the weight of current
        int start_weight = ComputeWeight(current,timeStamp);

        //copy of each transactions approvers
        std::vector<t_ptrTx> currentView = current->m_approvedBy;

        //store the weight of the sites that approuve current
        std::vector<int> sitesWeight;

        //store the probability to walk in each site that approuve current
        //pair = index,prob
        std::vector<std::pair<int,double>> sitesProb;

        //filter current View
        filterView(currentView,timeStamp);

        if(currentView.size() == 0)
        {
            break;
        }

        //if only one approver available dont compute the weight
        if(currentView.size() == 1)
        {
            current = currentView.front();
        }

        else
        {
            //choose if calculate the weights of available transaction here
            //random generator
            std::uniform_real_distribution<double> walkChoice(0.0,1.0);

            //the exp sum (see proba formula page 21)
            double sum_exp = 0.0;
            int weight;

            for(int j = 0; j < static_cast<int>(currentView.size()); j++)
            {
                weight = ComputeWeight(currentView[j],timeStamp);
                sum_exp = sum_exp + double(exp(double(-alphaVal*(start_weight - weight))));
                sitesWeight.push_back(weight);
            }

            //a probability (see formula page 21)
            double prob;

            for(int j = 0; j < static_cast<int>(sitesWeight.size()); j++)
            {
               prob = double(exp(double(-alphaVal*(start_weight - sitesWeight[j]))));
               prob = prob/sum_exp;
               sitesProb.push_back(std::make_pair(j,prob));
            }

            std::sort(sitesProb.begin(), sitesProb.end(),[](const std::pair<int,double> &a, const std::pair<int,double> &b){
            return a.second < b.second;}); //we sort sitesProb
            int nextCurrentIndex = 0; //index of the selected sites in CurrentView
            double probWalkChoice =  walkChoice(getTanglePtr()->getRandGen()); //the prob of the choosen site for the next walk

            for(int j = 0; j < static_cast<int>(sitesProb.size()) - 1; j++)
            {
                if(j == 0)
                {
                    if(probWalkChoice < sitesProb[j].second)
                    {
                        nextCurrentIndex = sitesProb[j].first;
                        break;
                    }
                }

                else
                {
                    if(probWalkChoice < sitesProb[j-1].second + sitesProb[j].second)
                    {
                        nextCurrentIndex = sitesProb[j].first;
                        break;
                    }
                }
            }

            current = currentView[nextCurrentIndex];
        }
    }

    walk_time = walkCounts;
    current->m_walkBacktracks = walkCounts;

    return current;
}

//initiate the TSA without the tracking
t_txApproved TxActor::WhitePaperTSA(double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int W, int N)
{
    //random generator with uniform distribution between W and 2*W for sites distance from tips
    std::uniform_int_distribution<int> backTrackDist(W,2*W);

    std::vector<t_ptrTx> start_sites; //vect of start sites for the TSA
    t_ptrTx temp_site;
    int backTD;

    //we pick randomly N sites for the TSA start point (may be the same)
    for(int i = 0; i < N; i++)
    {
        backTD = backTrackDist(getTanglePtr()->getRandGen());
        temp_site = getWalkStart(tips,backTrackDist(getTanglePtr()->getRandGen()));
        start_sites.push_back(temp_site);
    }

    t_txApproved selected_tips; //tips selected by the TSA (may be the same)
    std::vector<std::pair<int,int>> walk_total; //Execution time (i.e walk time) of the TSA for each start sites
    int walk_time; //Execution time for one start site

    //initiate for each start site a TSA and we save the execution time in time_exe
    for(int i = 0; i < N; i++)
    {
        if(i != N-1)
        {
            selected_tips.push_back(WPWalkTipSelection(start_sites[i],alphaVal,tips,timeStamp,walk_time));
            walk_total.push_back(std::make_pair(i,walk_time));
        }

        else
        {
            selected_tips.push_back(WPWalkTipSelection(start_sites[i],alphaVal,tips,timeStamp,walk_time));
            walk_total.push_back(std::make_pair(i,walk_time));
        }
    }

    std::sort(walk_total.begin(), walk_total.end(),[](const std::pair<int,int> &a, const std::pair<int,int> &b){
    return a.second < b.second;});

    t_txApproved tipstoApprove; //the two first tips selected by the TSA based on time execution
    //index in the time_exe vector of the first tip to be selected based on time execution
    int Index_first_tips = walk_total[0].first;
    tipstoApprove.push_back(selected_tips[Index_first_tips]); //we push it

   for(long unsigned int j = 1; j < walk_total.size(); j++)
   {
       Index_first_tips = walk_total[j].first;

       if(selected_tips[Index_first_tips]->TxNumber != tipstoApprove[0]->TxNumber)
       {
           tipstoApprove.push_back(selected_tips[Index_first_tips]);
           break;
       }
   }

    return tipstoApprove;
}

/**********************************************************G-IOTA*************************************************************************/

//perform the TSA based on G-IOTA
t_txApproved TxActor::GIOTA(double alphaVal, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int W, int N)
{
   auto chosenTips = WhitePaperTSA(alphaVal,tips,timeStamp,W,N); //we launch the weighted TSA to get two tips
   t_txApproved copyTips; //vector of tips

   for(auto it = tips.begin(); it != tips.end(); ++it) //we push left behind tips
   {
       auto tip = it->second;

       for(int j = 0; j < static_cast<int>(chosenTips.size()); j++)
       {
           if(chosenTips[j]->TxNumber != tip ->TxNumber)
           {
               copyTips.push_back(tip);
               break;
           }
       }
   }

   if(!copyTips.size()) //if the vector is empty after, no need to choose a third tips to approve
   {
       return chosenTips;
   }

   //we choose randomly a third tips to approve
   std::uniform_int_distribution<int> ThirdTipsIdx(0,static_cast<int>(copyTips.size()) - 1);
   auto idx = ThirdTipsIdx(getTanglePtr()->getRandGen());
   auto LeftTips = copyTips[idx];

   chosenTips.push_back(LeftTips);
   return chosenTips;
}

/**********************************************************E-IOTA*************************************************************************/

//perform the TSA based on E-IOTA
t_txApproved TxActor::EIOTA(double p1, double p2, std::map<int, t_ptrTx>& tips, omnetpp::simtime_t timeStamp, int W, int N)
{
    std::uniform_real_distribution<> rGen(0, 1);
    double lowAlpha = 0.1;
    double highAlpha = 5.0;
    auto r = rGen(getTanglePtr()->getRandGen()); //choose randomly r 
    t_txApproved chosenTips;

    if(r < p1)
    {
        chosenTips = WhitePaperTSA(0.0,tips,timeStamp,W,N);
    }

    else if(p1 <= r && r < p2)
    {
       chosenTips = WhitePaperTSA(lowAlpha,tips,timeStamp,W,N); 
    }

    else
    {
        chosenTips = WhitePaperTSA(highAlpha,tips,timeStamp,W,N);
    }

    return chosenTips;
}

