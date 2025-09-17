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
createInteraction(pFlow::cfModels::sup<true>, pFlow::stationaryGeometry);
createInteraction(pFlow::cfModels::sup<false>, pFlow::stationaryGeometry);

// rotationAxisMotionGeometry
createInteraction(pFlow::cfModels::sup<true>, pFlow::rotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::sup<false>, pFlow::rotationAxisMotionGeometry);

// vibratingMotionGeometry
createInteraction(pFlow::cfModels::sup<true>, pFlow::vibratingMotionGeometry);
createInteraction(pFlow::cfModels::sup<false>, pFlow::vibratingMotionGeometry);

// conveyorBeltMotionGeometry
createInteraction(pFlow::cfModels::sup<true>, pFlow::conveyorBeltMotionGeometry);
createInteraction(pFlow::cfModels::sup<false>, pFlow::conveyorBeltMotionGeometry);

// multiRotationAxisMotionGeometry
createInteraction(pFlow::cfModels::sup<true>, pFlow::multiRotationAxisMotionGeometry);
createInteraction(pFlow::cfModels::sup<false>, pFlow::multiRotationAxisMotionGeometry);	