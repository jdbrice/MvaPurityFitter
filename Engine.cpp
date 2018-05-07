

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PurityFitter.h"
#include "CutCurves.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	

	Logger::setGlobalLogLevel( "none" );

	loguru::init(argc, argv);
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);
	
	TaskFactory::registerTaskRunner<PurityFitter>( "PurityFitter" );
	TaskFactory::registerTaskRunner<CutCurves>( "CutCurves" );
	

	TaskEngine engine( argc, argv, "PurityFitter" );


	return 0;
}
