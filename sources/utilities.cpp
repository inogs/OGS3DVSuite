#include "OGS/utilities.h"

#include <vtkFieldData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>

#include <chrono>

void OGS::utilities::strsplit(const std::string &str,
                              const std::string &splitBy,
                              std::vector<std::string> &tokens) {
  // Store the original string in the array, so we can loop the rest of the
  // algorithm.
  tokens.push_back(str);

  // Store the split index in a 'size_t' (unsigned integer) type.
  size_t splitAt;
  // Store the size of what we're splicing out.
  size_t splitLen = splitBy.size();
  // Create a string for temporarily storing the fragment we're processing.
  std::string frag;
  // Loop infinitely - break is internal.
  while (true) {
    // Store the last string in the vector, which is the only logical candidate
    // for processing.
    frag = tokens.back();
    // The index where the split is.
    splitAt = frag.find(splitBy);
    // If we didn't find a new split point...
    if (splitAt == std::string::npos)
      break;  // Break the loop and (implicitly) return.
    // Put everything from the left side of the split where the string being
    // processed used to be.
    tokens.back() = frag.substr(0, splitAt);
    // Push everything from the right side of the split to the next empty index
    // in the vector.
    tokens.push_back(
        frag.substr(splitAt + splitLen, frag.size() - (splitAt + splitLen)));
  }
}


void OGS::utilities::BuildTimeList(OGS::Time::TimeList &TL,
                                   vtkInformation *Info) {
  // Build a TimeList object using the pipeline temporal data.
  // This TimeList will be later used for computing the averages
  // given a TimeRequestor. The metadata array might not be available
  // from the beginning.

  int ntsteps = Info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  double *tsteps;
  tsteps = Info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  // Create a TimeObjectList where to store all the instants
  OGS::Time::TimeObjectList TOL(ntsteps);

  // Iterate the number of steps and set the values of the list
  for (int ii = 0; ii < ntsteps; ++ii) {
    // Convert to struct tm
    time_t time = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::time_point(
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::duration<double>(tsteps[ii]))));
    struct tm tm = *localtime(&time);

    // Set up the TimeObjectList
    TOL[ii] = OGS::Time::TimeObject(tm);
  }

  // Sort and print the list (for debugging purposes)
  TOL.sort();
  //	printf("TOL: %s\n",TOL.as_string("%Y%m%d").c_str());

  // Now create the TimeList from the TimeObjectList
  TL = OGS::Time::TimeList(TOL);
  //	printf("Defined list of %d elements: %s ...
  //%s\n",TL.len(),TL[0].as_string("%Y-%m-%d
  //%H:%M:%S").c_str(),TL[-1].as_string("%Y-%m-%d %H:%M:%S").c_str());
}


void OGS::utilities::RecoverMasterFileName(std::string &fname,
                                           vtkDataSet *input) {
  // Recover the master file name from the metadata array
  // Return whether we need to stop executing or not

  vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
      input->GetFieldData()->GetAbstractArray("Metadata"));

  // If successful, recover the file name
  if (vtkmetadata)
    fname = vtkmetadata->GetValue(4);
  else
    fname = std::string("NotFound");
}
