

#include "itkFastGrowCut.h"

#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include "itkTestingMacros.h"
// class ShowProgress : public itk::Command
// {
// public:
//   itkNewMacro(ShowProgress);

//   void
//   Execute(itk::Object * caller, const itk::EventObject & event) override
//   {
//     Execute((const itk::Object *)caller, event);
//   }

//   void
//   Execute(const itk::Object * caller, const itk::EventObject & event) override
//   {
//     if (!itk::ProgressEvent().CheckEvent(&event))
//     {
//       return;
//     }
//     const auto * processObject = dynamic_cast<const itk::ProcessObject *>(caller);
//     if (!processObject)
//     {
//       return;
//     }
//     std::cout << " " << processObject->GetProgress();
//   }
// };

int
main(int argc, char * argv[])
{
  const char * our_initial = "../out/initial.mha";
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
    ImageType::IndexType corner = {45, 45, 45};
    ImageType::SizeType cubeSize = {10, 10, 10};
    ImageType::RegionType cubeRegion(corner, cubeSize);
    itk::ImageRegionIterator<ImageType> it(image, cubeRegion);
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        std::cout << "Position (x, y, z): " << it.GetIndex() << std::endl;
        it.Set(1);
    }

    ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(image, our_initial, true));

    // seed at one corner of the volume for one class and the center of the volume for the other class
    LabelType::Pointer seeds = LabelType::New();
    seeds->SetRegions(region);
    seeds->Allocate();
    seeds->FillBuffer(0);
    LabelType::IndexType cornerSeed = {0, 0, 0};
    seeds->SetPixel(cornerSeed, 1);
    LabelType::IndexType centerSeed = {50, 50, 50};
    seeds->SetPixel(centerSeed, 20);

  using FGCType = itk::FastGrowCut<ImageType, LabelType>;
  FGCType::Pointer fgcFilter = FGCType::New();

  ITK_EXERCISE_BASIC_OBJECT_METHODS(fgcFilter, FastGrowCut, ImageToImageFilter);

  // Read the images
//   /ImageType::Pointer image = itk::ReadImage<ImageType>(inputImage);
  //LabelType::Pointer seeds = itk::ReadImage<LabelType>(seedsImage);

//   ShowProgress::Pointer showProgress = ShowProgress::New();
//   fgcFilter->AddObserver(itk::ProgressEvent(), showProgress);

  // Filter the original, possibly noisy image
  fgcFilter->SetInput(image);
  fgcFilter->SetSeedImage(seeds);
  fgcFilter->Update();
  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(fgcFilter->GetOutput(), out_processed, true));

  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
