/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment,
** also known as GENESIS 3 base code.
**           copyright (C) 2003-2009 Upinder S. Bhalla. and NCBS
** It is made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#ifndef _StateScanner_h
#define _StateScanner_h
class StateScanner
{
#ifdef DO_UNIT_TESTS
	friend void testStateScanner();
#endif
	public:
		StateScanner();
		~StateScanner();
		
		///////////////////////////////////////////////////
		// Field function definitions
		///////////////////////////////////////////////////
		static double getSettleTime( Eref e );
		static void setSettleTime( const Conn* c, double value );
		static double getSolutionSeparation( Eref e );
		static void setSolutionSeparation( const Conn* c, double value );
		static unsigned int getStateCategory(Eref e, const unsigned int& i);
		unsigned int localGetStateCategory( unsigned int i) const;
		static void setStateCategory( 
			const Conn* c, unsigned int val, const unsigned int& i );
		static void addTrackedMolecule( const Conn* c, Id val );
		static void dropTrackedMolecule( const Conn* c, Id val );
		/*
		static unsigned int getNumTrackedMolecules( Eref e );
		static void setNumTrackedMolecules( const Conn* c, unsigned int value );
		void localSetNumTrackedMolecules( unsigned int value );
		static Id getTrackedMolecule( Eref e, const unsigned int& i );
		static void setTrackedMolecule(
			const Conn* c, Id val, const unsigned int& i );
		Id localGetTrackedMolecule( unsigned int i ) const;
		void localSetTrackedMolecule( Id elm, unsigned int i );
		*/
		static unsigned int getClassification( Eref e );

		///////////////////////////////////////////////////
		// Msg Dest function definitions
		///////////////////////////////////////////////////
		
		static void doseResponse( const Conn* c, 
			Id variableMol, 
			double start, double end, 
			unsigned int numSteps );

		static void logDoseResponse( const Conn* c, 
			Id variableMol, 
			double start, double end, 
			unsigned int numSteps );

		void innerDoseResponse( Id variableMol, 
			double start, double end, 
			unsigned int numSteps,
			bool useLog );

		static void classifyStates( const Conn* c, 
			unsigned int numStartingPoints,
			bool useMonteCarlo,
			bool useLog );
		void innerClassifyStates(
			unsigned int numStartingPoints,
			bool useMonteCarlo,
			bool useLog );

		// funcs to handle externally imposed changes in mol N
		static void setMolN( const Conn* c, double y, unsigned int i );
		static void assignStoichFunc( const Conn* c, void* stoich );
		void assignStoichFuncLocal( void* stoich );

		
	private:
		bool isMoleculeIndexGood( unsigned int i ) const;
		///////////////////////////////////////////////////
		// Internal fields.
		///////////////////////////////////////////////////
		double settleTime_;
		double solutionSeparation_;
		unsigned int numTrackedMolecules_;
		vector< Id > trackedMolecule_;
		vector< unsigned int> stateCategories_;

		unsigned int numSolutions_;
		unsigned int numStable_;
		unsigned int numSaddle_;
		unsigned int numOsc_;
		unsigned int numOther_;
		unsigned int classification_;
};

extern const Cinfo* initStateScannerCinfo();
#endif // _StateScanner_h
