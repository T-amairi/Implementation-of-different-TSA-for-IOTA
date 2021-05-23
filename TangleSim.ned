//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//

simple TxActorModule
{
    parameters:
        volatile double txGenRate @unit(s) = 3.0s; // how often a transactor will issue a transaction
        volatile double powTime @unit(s) = 10.0s; // time taken to compute proof of work to approve two transactions

        //determines randomness of tip selection walk, higher is more deterministic, lower is more random
        double walkAlphaValue = 0.8; // *******MUST BE BETWEEN 0 AND 1********* 
        int walkDepth = 5;

        //determines which tip selection + attach method to use ()
        string tipSelectionMethod = "WPWALKTrack";

        //WPWALK && GIOTA && EIOTA
        double Alpha = 5.0;
        int WalkDistance = 250;
        int NumberofWalkers = 100;

        //EIOTA
        double p1 = 0.1;
        double p2 = 0.65;

        bool recordWeights = default(false);

        // Denotes how many walkers we use in WALK, k being the number of tips selected per attach
        int k_Multiplier = 2;

    gates:
        inout tangleConnect;

}

simple TangleModule
{
    parameters:
        @display("i=block/routing" );
        int transactionLimit = 100; // how many transactions to simulate before stopping

        string tipDataFilename = default(" "); //default("./data/RTracking/GeneraTipData.txt");
        string tipAgeFilename = default(" "); //default("./data/RTracking/TipAge.txt");
        string blockWeightFilename = default(" ");//default("./data/RTracking/BlockWeight.txt");

    gates:
        inout actorConnect[];

}

network TangleSim
{
    parameters:
        int txActorNumber = 2; // How many transactors in our simulation - must be assigned in .ini file

        @display("bgb=376,219");
    submodules:
        actors[txActorNumber]: TxActorModule {
            @display("p=110,142");
        }
        tangle: TangleModule {
            @display("p=314,142");
        }
    connections:
        for i=0..txActorNumber -1 {
            actors[i].tangleConnect <--> {  delay = 1ms; } <--> tangle.actorConnect++;
        }

}
