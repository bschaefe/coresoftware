#include "ActsGeometry.h"
#include "alignmentTransformationContainer.h"
#include "sPHENIXActsDetectorElement.h"


sPHENIXActsDetectorElement::~sPHENIXActsDetectorElement() = default;

const Acts::Transform3& sPHENIXActsDetectorElement::transform(const Acts::GeometryContext& ctxt)  const
{

  if(alignmentTransformationContainer::use_alignment)
    {
      Acts::GeometryIdentifier id = surface().geometryId();

      unsigned int volume = id.volume(); 
      unsigned int layer = id.layer(); 
      unsigned int sphlayer = base_layer_map.find(volume)->second + layer / 2 -1;
      unsigned int sensor = id.sensitive() - 1;  // Acts sensor ID starts at 1

      const std::vector<std::vector<Acts::Transform3>>& transformVec 
	= ctxt.get<std::vector<std::vector<Acts::Transform3>>&>();
            
      auto& layerVec = transformVec[sphlayer];    // get the vector of transforms for this layer
      if(layerVec.size() > sensor)
	{
	  //std::cout << "sPHENIXActsDetectorElement: return transform:  volume " << volume <<" Acts  layer " << layer << " sensor " << sensor
	  //	  	    << " sphenix layer " << sphlayer << " layerVec size " << layerVec.size() << std::endl;
	
	  return layerVec[sensor];
	}
      
      // if we are still here, it was not found
      std::cout << " Alignment transform not found, for identifier " << id << " use construction transform " << std::endl;
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      //std::cout << "           construction transform: " << std::endl << transform.matrix() << std::endl;      
      return transform;
    }
 else
    {
      // return the construction transform
      const Acts::Transform3& transform = TGeoDetectorElement::transform(ctxt);  // ctxt is unused here
      return transform;
    }

}


