#include "sphereInteraction.hpp"
#include "geometryMotions.hpp"
#include "contactForceModels.hpp"
#include "unsortedContactList.hpp"
#include "sortedContactList.hpp"
#include "createBoundarySphereInteraction.hpp"



#define createInteraction(ForceModel,GeomModel) 	\
													\
	template class pFlow::sphereInteraction< 		\
		ForceModel,									\
		GeomModel,									\
		pFlow::unsortedContactList>;				\
													\
	template class pFlow::sphereInteraction< 		\
		ForceModel,									\
		GeomModel,									\
		pFlow::sortedContactList>;					\
	createBoundarySphereInteraction(ForceModel, GeomModel)

// ============ SUP 模型 ============

// stationaryGeometry
createInteraction(pFlow::cfModels::limitedSUPRollingCDT, pFlow::stationaryGeometry);
createInteraction(pFlow::cfModels::nonLimitedSUPRollingCDT, pFlow::stationaryGeometry);

// rotationAxisMotionGeometry
createInteraction(pFlow::cfModels::limitedSUPRollingCDT, pFlow::rotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedSUPRollingCDT, pFlow::rotationAxisMotionGeometry);

// vibratingMotionGeometry
createInteraction(pFlow::cfModels::limitedSUPRollingCDT, pFlow::vibratingMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedSUPRollingCDT, pFlow::vibratingMotionGeometry);

// conveyorBeltMotionGeometry
createInteraction(pFlow::cfModels::limitedSUPRollingCDT, pFlow::conveyorBeltMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedSUPRollingCDT, pFlow::conveyorBeltMotionGeometry);

// multiRotationAxisMotionGeometry
createInteraction(pFlow::cfModels::limitedSUPRollingCDT, pFlow::multiRotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::nonLimitedSUPRollingCDT, pFlow::multiRotationAxisMotionGeometry);