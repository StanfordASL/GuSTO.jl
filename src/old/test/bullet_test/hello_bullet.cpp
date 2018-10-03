#include <iostream>
#include <sstream>
#include <fstream>

#include "btBulletCollisionCommon.h"
#include "BulletCollision/CollisionShapes/btConvex2dShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkPairDetector.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpaPenetrationDepthSolver.h"
#include "BulletCollision/NarrowPhaseCollision/btPointCollector.h"
#include "BulletDynamics/Featherstone/btMultiBody.h"
#include "LinearMath/btConvexHull.h"
#include <btBulletDynamicsCommon.h>

int main() {
  btCollisionConfiguration* m_coll_config   = new btDefaultCollisionConfiguration();
  btCollisionDispatcher* m_dispatcher       = new btCollisionDispatcher(m_coll_config);
  btBroadphaseInterface* m_broadphase       = new btDbvtBroadphase();
  btCollisionWorld* m_world                 = new btCollisionWorld(m_dispatcher, m_broadphase, m_coll_config);

  btCollisionObject* box = new btCollisionObject();
  btBoxShape* box_shape = new btBoxShape(btVector3(5.0,5.0,5.0));
  box_shape->setMargin(0.);
  box->getWorldTransform().setOrigin(btVector3(0.0,0.0,0.0));
  box->setCollisionShape(box_shape);
  
  m_world->addCollisionObject(box);
  m_world->updateAabbs();

  m_world->getNumCollisionObjects();
  delete m_world;
  delete m_broadphase;
	delete m_dispatcher;
	delete m_coll_config;
  
  printf("Out safely \n");
  return 0;
}
