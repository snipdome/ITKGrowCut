

#include "itkFastGrowCut.h"
#include "itkFastGrowCut.cxx"

#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include "itkTestingMacros.h"

namespace
{
class ShowProgress : public itk::Command
{
public:
  itkNewMacro(ShowProgress);

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * caller, const itk::EventObject & event) override
  {
    if (!itk::ProgressEvent().CheckEvent(&event))
    {
      return;
    }
    const auto * processObject = dynamic_cast<const itk::ProcessObject *>(caller);
    if (!processObject)
    {
      return;
    }
    std::cout << " " << processObject->GetProgress();
  }
};
}

int
main(int argc, char * argv[])
{
  const char * our_initial = "../out/initial.mha";
  const char * our_seed = "../out/seed.mha";
  const char * out_processed = "../out/processed.mha";

  constexpr unsigned int Dimension = 3;
  using PixelType = short;
  using ImageType = itk::Image<PixelType, Dimension>;
  using LabelType = itk::Image<unsigned char, Dimension>;

  // inputImage is a 100x100x100 with a small cube in the middle
  // fill with 0, but fill the middle cube with 1
  ImageType::Pointer image = ImageType::New();
  ImageType::RegionType region;
  ImageType::IndexType start;
  start.Fill(0);
  ImageType::SizeType size;
  size.Fill(100);
  region.SetSize(size);
  region.SetIndex(start);
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);
  ImageType::IndexType corner = {30, 30, 30};
  ImageType::SizeType cubeSize = {40, 40, 40};
  ImageType::RegionType cubeRegion(corner, cubeSize);
  itk::ImageRegionIterator<ImageType> it(image, cubeRegion);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) it.Set(100);
  itk::WriteImage(image, our_initial, true);


  // Seed image is a 100x100x100 with a seed at one corner of the volume for one class and the center of the volume for the other class
  LabelType::Pointer seed_image = LabelType::New();
  seed_image->SetRegions(region);
  seed_image->Allocate();
  seed_image->FillBuffer(0);
  LabelType::IndexType cornerSeed = {0, 0, 0};
  seed_image->SetPixel(cornerSeed, 1);
  // LabelType::IndexType centerSeed = {50, 50, 50};
  // seed_image->SetPixel(centerSeed, 2);
  itk::WriteImage(seed_image,  our_seed, true);


  using FGCType = itk::FastGrowCut<ImageType, LabelType>;
  FGCType::Pointer fgcFilter = FGCType::New();

  ShowProgress::Pointer showProgress = ShowProgress::New();
  fgcFilter->AddObserver(itk::ProgressEvent(), showProgress);

  // Filter the original, possibly noisy image
  fgcFilter->SetInput(image);
  fgcFilter->SetSeedImage(seed_image);
  fgcFilter->Update();
  itk::WriteImage(fgcFilter->GetOutput(), out_processed, true);

  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
