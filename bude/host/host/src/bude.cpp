#include "CL/opencl.h"
//#include "CL/opencl.h"
//#include "CL/opencl.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
//#include <Winsock2.h>
//#include "AOCL_Utils.h"
#include "AOCLUtils/aocl_utils.h"



using namespace aocl_utils;

#define POSES_PER_WI  4
#define WGSIZE       256


#define ALTERA_OPENCL


// OpenCL runtime configuration
cl_platform_id platform = NULL;
//unsigned num_devices = 0;
cl_device_id device; // num_devices elements
cl_context context = NULL;
cl_command_queue queue; // num_devices elements
cl_program program = NULL;
cl_kernel kernel; // num_devices elements


/*
struct
{
  cl_platform_id   platform;
  cl_device_id     device;
  cl_context       context;
  cl_command_queue queue;
  cl_program       program;
  cl_kernel        kernel;
} cl;
*/

bool initCL();
void releaseCL();
unsigned num_devices = 1;
cl_uint num_platforms;


typedef struct
{
  cl_float x, y, z;
  cl_int type;
} Atom;

typedef struct
{
  cl_int   hbtype;
  cl_float radius;
  cl_float hphb;
  cl_float elsc;
} FFParams;


#include "kernel.h"
#include "protein.h"
#include "ligand.h"
#include "forcefield.h"
#include "poses.h"

static void dump_error(const char *str, cl_int status) {
	printf("%s\n", str);
	printf("Error code: %d\n", status);
}

// free the resources allocated during initialization
static void freeResources() {

/*	if (kernel)
		clReleaseKernel(kernel);
	if (kernel_read)
		clReleaseKernel(kernel_read);
	if (kernel_write)
		clReleaseKernel(kernel_write);
	if (program)
		clReleaseProgram(program);
	if (queue)
		clReleaseCommandQueue(queue);
	if (hdatain)
		clSVMFreeAltera(context, hdatain);
	if (hdataout)
		clSVMFreeAltera(context, hdataout);
	if (context)
		clReleaseContext(context);*/

}



int main(int argc, char *argv[])
{
  printf("Starting BUDE Altera\n");
  cl_int err;

  //memset(&cl, 0, sizeof(cl));
  initCL();

  cl_int natlig = sizeof(h_ligand)/sizeof(Atom);
  cl_int natpro = sizeof(h_protein)/sizeof(Atom);
  cl_int ntypes = sizeof(h_forcefield)/sizeof(FFParams);
  //cl_int nposes = sizeof(h_poses0)/sizeof(float);
  //cl_int nposes = sizeof(h_poses0)/sizeof(float);
  cl_int nposes = 256; 
  cl_mem protein, ligand, energies, forcefield, poses[6];

  printf("nposes are %d\n",nposes);

  printf("Creating buffers\n");
  // Create buffers

  protein = (cl_mem)clSVMAllocAltera(context, 0, natpro*sizeof(Atom), 1024);
  if (!protein) {
	//  dump_error("Failed to allocate buffers.", status);
	  //freeResources();
	  return 1;

  }

  /*protein = clCreateBuffer(context, CL_MEM_READ_ONLY,
                           natpro*sizeof(Atom), NULL, &err);
  checkError(err, "creating protein");*/


  ligand = (cl_mem)clSVMAllocAltera(context, 0, natlig*sizeof(Atom), 1024);
  if (!ligand) {
	  //  dump_error("Failed to allocate buffers.", status);
	  //freeResources();
	  return 1;

  }

  /*ligand = clCreateBuffer(context, CL_MEM_READ_ONLY,
                          natlig*sizeof(Atom), NULL, &err);
  checkError(err, "creating ligand");*/

  energies = (cl_mem)clSVMAllocAltera(context, 0, nposes*sizeof(cl_float), 1024);
  if (!energies) {
	  //  dump_error("Failed to allocate buffers.", status);
	  //freeResources();
	  return 1;

  }

  /*
  energies = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                            nposes*sizeof(cl_float), NULL, &err);
  checkError(err, "creating energies");

  */

  forcefield = (cl_mem)clSVMAllocAltera(context, 0, ntypes*sizeof(FFParams), 1024);
  if (!forcefield) {
	  //  dump_error("Failed to allocate buffers.", status);
	  //freeResources();
	  return 1;
  }
  /*
  forcefield = clCreateBuffer(context, CL_MEM_READ_ONLY,
                              ntypes*sizeof(FFParams), NULL, &err);
  checkError(err, "creating forcefield");
  */

 /* global_forcefield = (cl_mem)clSVMAllocAltera(context, 0, ntypes*sizeof(FFParams), 1024);
  if (!global_forcefield) {
	  //  dump_error("Failed to allocate buffers.", status);
	  //freeResources();
	  return 1;
  }*/

  /*global_forcefield = clCreateBuffer(context, CL_MEM_READ_ONLY,
	  ntypes*sizeof(FFParams), NULL, &err);
  checkError(err, "creating global forcefield");*/
  


  for (int i = 0; i < 6; i++)
  {
	  poses[i] = (cl_mem)clSVMAllocAltera(context, 0, nposes*sizeof(cl_float), 1024);
	  if (!poses[i]) {
		  //  dump_error("Failed to allocate buffers.", status);
		  //freeResources();
		  return 1;
	  }
	  /*poses[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                              nposes*sizeof(cl_float), NULL, &err);*/
  }
  
  printf("SVM memory allocated. No writting data to device\n");
  // Write data to device

  /*err = clEnqueueWriteBuffer(queue, protein, CL_TRUE,
        0,natpro*sizeof(Atom), h_protein, 0, NULL, NULL);

  checkError(err, "writing protein");

  err = clEnqueueWriteBuffer(queue, ligand, CL_TRUE, 0,
                             natlig*sizeof(Atom), h_ligand,
                             0, NULL, NULL);
  checkError(err, "writing ligand");
  err = clEnqueueWriteBuffer(queue, forcefield, CL_TRUE, 0,
                             ntypes*sizeof(FFParams), h_forcefield,
                             0, NULL, NULL);
  checkError(err, "writing forcefield");

  err = clEnqueueWriteBuffer(queue, poses[0], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses0,
                             0, NULL, NULL);
  checkError(err, "writing poses 0");
  err = clEnqueueWriteBuffer(queue, poses[1], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses1,
                             0, NULL, NULL);
  checkError(err, "writing poses 1");
  err = clEnqueueWriteBuffer(queue, poses[2], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses2,
                             0, NULL, NULL);
  checkError(err, "writing poses 2");
  err = clEnqueueWriteBuffer(queue, poses[3], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses3,
                             0, NULL, NULL);
  checkError(err, "writing poses 3");
  err = clEnqueueWriteBuffer(queue, poses[4], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses4,
                             0, NULL, NULL);
  checkError(err, "writing poses 4");
  err = clEnqueueWriteBuffer(queue, poses[5], CL_TRUE, 0,
                             nposes*sizeof(cl_float), h_poses5,
                             0, NULL, NULL);
  checkError(err, "writing poses 5");
  */
  
  printf("Set kernel arguments\n");

  // Set kernel arguments
  //err = clSetKernelArgSVMPointerAltera(kernel, 0, (void*)natlig);
  err  = clSetKernelArg(kernel, 0, sizeof(cl_int), &natlig);
  //err |= clSetKernelArgSVMPointerAltera(kernel, 1, (void*)natpro);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_int), &natpro);
  err |= clSetKernelArgSVMPointerAltera(kernel, 2, (void*)protein);
  //err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &protein);
  err |= clSetKernelArgSVMPointerAltera(kernel, 3, (void*)ligand);
  //err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &ligand);
  err |= clSetKernelArgSVMPointerAltera(kernel, 4, (void*)(poses[0]));
  //err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), poses+0);
  err |= clSetKernelArgSVMPointerAltera(kernel, 5, (void*)(poses[1]));
  //err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), poses+1);
  err |= clSetKernelArgSVMPointerAltera(kernel, 6, (void*)(poses[2]));
  //err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), poses+2);
  err |= clSetKernelArgSVMPointerAltera(kernel, 7, (void*)(poses[3]));
  //err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), poses+3);
  err |= clSetKernelArgSVMPointerAltera(kernel, 8, (void*)(poses[4]));
  //err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), poses+4);
  err |= clSetKernelArgSVMPointerAltera(kernel, 9, (void*)(poses[5]));
  //err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), poses+5);
  err |= clSetKernelArgSVMPointerAltera(kernel, 10, (void*)(energies));
  //err |= clSetKernelArg(kernel, 10, sizeof(cl_mem), &energies);
  err |= clSetKernelArgSVMPointerAltera(kernel, 11, (void*)forcefield);
  //err |= clSetKernelArg(kernel, 11, sizeof(cl_mem), &forcefield);
  //err |= clSetKernelArgSVMPointerAltera(kernel, 12, (void*)global_forcefield);
  //err |= clSetKernelArg(kernel, 12, ntypes*sizeof(FFParams), NULL);
  //err |= clSetKernelArg(kernel, 12, ntypes*sizeof(FFParams), global_forcefield);
  //err |= clSetKernelArgSVMPointerAltera(kernel, 13, (void*)ntypes);
  err |= clSetKernelArg(kernel, 12, sizeof(cl_int), &ntypes);
  //err |= clSetKernelArgSVMPointerAltera(kernel, 14, (void*)nposes);
  err |= clSetKernelArg(kernel, 13, sizeof(cl_int), &nposes);
  checkError(err, "setting arguments");

  size_t global[1] = {nposes/POSES_PER_WI};
  size_t local[1] = {WGSIZE};

  cl_event pro_event;
  
  printf("Doing warming-up\n");

//  double start = getCurrentTimestamp() * 1.0e+9; //ns
  // Warm-up run (not timed)

  err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
	  global, local, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
	  dump_error("Failed clEnqueueNDRangeKernel", err);
	  freeResources();
	  return 1;
  }

  //checkError(err, "queuing kernel");
  //printf("Running kernel\n");
  err = clFinish(queue);
  checkError(err, "running kernel");
//  double end = getCurrentTimestamp() * 1.0e+9; //ns

//  double ms = ((end-start))*1.0e-6;
//  printf("Warm up time = %.1lf ms\n", ms);
  
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double start = tv.tv_usec + tv.tv_sec*1e6;

  /*SYSTEMTIME      tv;
  GetLocalTime( &tv );   
  printf( "Started at %02d:%02d:%02d:%06d\n", tv.wHour,    tv.wMinute, tv.wSecond, tv.wMilliseconds * 1000 );
  double start = tv.wMinute*60*1000 + tv.wSecond*1000 + tv.wMilliseconds;*/

//  start = getCurrentTimestamp() * 1.0e+9; //ns
  
  
  printf("Doing timed runs\n");

  // Timed runs
#define RUNS 10
  for (int i = 0; i < RUNS; i++)
  {
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
                                 global, local, 0, NULL, &pro_event);
	//clGetProfileInfoAltera(pro_event);
  }

  err = clFinish(queue);
  checkError(err, "running kernel");


  gettimeofday(&tv, NULL);
  double end = tv.tv_usec + tv.tv_sec*1e6;
  //end = getCurrentTimestamp() * 1.0e+9; //ns
  
  double ms = ((end-start)/RUNS)*1.0e-6;
  printf("Average run time = %.1lf ms\n", ms);

  /*GetLocalTime( &tv );   
  printf( "End at %02d:%02d:%02d:%06d\n", tv.wHour,    tv.wMinute, tv.wSecond, tv.wMilliseconds * 1000 );
  double end = tv.wMinute*60*1000 + tv.wSecond*1000 + tv.wMilliseconds;*/


  printf("Reading results\n");

  // Read results
  float *results = (float *)(malloc(nposes*sizeof(float)));
  err = clEnqueueReadBuffer(queue, energies, CL_TRUE, 0,
                            nposes*sizeof(cl_float), results,
                            0, NULL, NULL);
  checkError(err, "reading results");
  for (int i = 0; i < 8; i++)
  {
    printf("Energy %d = %.2f\n", i, results[i]);
  }
  free(results);

  ms = ((end-start)/RUNS)*1.0e-6;
  printf("Average time = %.1lf ms\n", ms);

  // Compute FLOP/s
  double runtime   = ms*1e-3;
  double poses_per_wi = POSES_PER_WI;
  double ops_per_wi = 27*poses_per_wi
    + natlig*(3 + 18*poses_per_wi + natpro*(11 + 30*poses_per_wi))
    + poses_per_wi;
  double total_ops     = ops_per_wi * (nposes/poses_per_wi);
  double flops      = total_ops / runtime;
  double gflops     = flops / 1e9;

  double interactions         = nposes * natlig * natpro;
  double interactions_per_sec = interactions / runtime;
  double work_items_per_sec = (nposes/poses_per_wi)/runtime;

  // Print stats
  printf("work_items/s:  %5.3lf million\n", (work_items_per_sec / 1e6));
  printf("Interactions/s:  %5.2lf billion\n", (interactions_per_sec / 1e9));
  printf("GFLOP/s:         %6.3lf\n", gflops);

  clReleaseMemObject(protein);
  clReleaseMemObject(ligand);
  clReleaseMemObject(energies);
  clReleaseMemObject(forcefield);
  clReleaseMemObject(poses[0]);
  clReleaseMemObject(poses[1]);
  clReleaseMemObject(poses[2]);
  clReleaseMemObject(poses[3]);
  clReleaseMemObject(poses[4]);
  clReleaseMemObject(poses[5]);

  releaseCL();
}

#ifdef ALTERA_OPENCL

bool initCL()
{
  
  cl_int status;

  printf("Initializing OpenCL\n");

  if(!setCwdToExeDir()) {
    return false;
  }

  // get the platform ID
  status = clGetPlatformIDs(1, &platform, &num_platforms);
  if (status != CL_SUCCESS) {
	  dump_error("Failed clGetPlatformIDs.", status);
	  freeResources();
	  return 1;
  }
  if (num_platforms != 1) {
	  printf("Found %d platforms!\n", num_platforms);
	  freeResources();
	  return 1;
  }

  // get the device ID
  status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, &num_devices);
  if (status != CL_SUCCESS) {
	  dump_error("Failed clGetDeviceIDs.", status);
	  freeResources();
	  return 1;
  }
  if (num_devices != 1) {
	  printf("Found %d devices!\n", num_devices);
	  freeResources();
	  return 1;
  }

  // create a context
  context = clCreateContext(0, 1, &device, NULL, NULL, &status);
  if (status != CL_SUCCESS) {
	  dump_error("Failed clCreateContext.", status);
	  freeResources();
	  return 1;
  }

  // Get the OpenCL platform.
/*  platform = findPlatform("Intel");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel FPGA OpenCL platform.\n");
    return false;
  }

 // Query the available OpenCL devices.
  scoped_array<cl_device_id> devices;
  cl_uint num_devices;

  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

  // We'll just use the first device.
  device = devices[0];

  printf("Platform: %s\n", getPlatformName(platform).c_str());
  printf("Using %d device(s)\n", num_devices);
  if (num_devices > 1)
  {
	printf("Error : Too many devices\n");
	exit(1);
  }
   
  //for(unsigned i = 0; i < num_devices; ++i) {
  printf("  %s\n", getDeviceName(device).c_str());
  //}

  // Create the context.
  context = clCreateContext(NULL, num_devices, &device, NULL, NULL, &status);
  checkError(status, "Failed to create context");
  */
  printf("Context created\n");

  /*
  size_t binsize = 0;
  unsigned char * binary_file = loadBinaryFile("bin/mem_bandwidth.aocx", &binsize);

  if (!binary_file) {
	  dump_error("Failed loadBinaryFile.", status);
	  freeResources();
	  return 1;
  }
  program = clCreateProgramWithBinary(context, 1, &device, &binsize, (const unsigned char**)&binary_file, &kernel_status, &status);
  
 */

  // Create the program for all device. Use the first device as the
  // representative device (assuming all device are of the same type).
  std::string binary_file = getBoardBinaryFile("budeMultiTD", device);
  printf("Using AOCX: %s\n", binary_file.c_str());
  program = createProgramFromBinary(context, binary_file.c_str(), &device, num_devices);

  // Build the program that was just created.
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  // Create per-device objects.
  //queue.reset(num_devices);
  //kernel.reset(num_devices);


  //for(unsigned i = 0; i < num_devices; ++i) {
    // Command queue.
    queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
    checkError(status, "Failed to create command queue");

    // Kernel.
    const char *kernel_name = "fasten_main";
    kernel = clCreateKernel(program, kernel_name, &status);
    checkError(status, "Failed to create kernel");

  //}

  printf("Initialization completed\n");

  return true;
}

#else

void initCL()
{

#define PLATFORM 0
#define DEVICE   0
#define MAX_NUM  8

  cl_int err;

  cl_platform_id platform;
  cl_platform_id platforms[MAX_NUM];
  err = clGetPlatformIDs(MAX_NUM, platforms, NULL);
  checkError(err, "getting platforms");
  platform = platforms[PLATFORM];

  cl_device_id devices[MAX_NUM];
  err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, MAX_NUM, devices, NULL);
  checkError(err, "getting devices");
  cl.device = devices[DEVICE];

  char name[128];
  clGetDeviceInfo(device, CL_DEVICE_NAME, 128, name, 0);
  printf("Using device: %s\n", name);

  cl.context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
  checkError(err, "creating context");

  queue = clCreateCommandQueue(
    context, device, CL_QUEUE_PROFILING_ENABLE, &err);
  checkError(err, "creating queue");

  program = clCreateProgramWithSource(
    context, 1, (const char**)&bude_kernel, NULL, &err);
  checkError(err, "creating program");

  char options[256];
  sprintf(options,
          "-cl-fast-relaxed-math -cl-mad-enable -DNUM_TD_PER_THREAD=%d",
          POSES_PER_WI);
  err = clBuildProgram(cl.program, 1, &cl.device, options, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    if (err == CL_BUILD_PROGRAM_FAILURE)
    {
      char log[16384];
      clGetProgramBuildInfo(cl.program, cl.device, CL_PROGRAM_BUILD_LOG,
                            16384, log, NULL);
      printf("%s\n", log);
    }
  }
  checkError(err, "building program");

  cl.kernel = clCreateKernel(cl.program, "fasten_main", &err);
  checkError(err, "creating kernel");
}

#endif

#define RELEASE(func, obj) if (obj) {func(obj); obj=NULL;};
void releaseCL()
{
  RELEASE(clReleaseKernel, kernel);
  RELEASE(clReleaseProgram, program);
  RELEASE(clReleaseCommandQueue, queue);
  RELEASE(clReleaseContext, context);
}

#ifdef ALTERA_OPENCL
#else

void checkError(cl_int err, const char *op)
{
  if (err != CL_SUCCESS)
  {
    printf("Error during operation '%s' (%d)\n", op, err);
    releaseCL();
  }
}
#endif

#ifdef ALTERA_OPENCL
// Free the resources allocated during initialization
void cleanup() {
  for(unsigned i = 0; i < num_devices; ++i) {
    if(kernel) {
      clReleaseKernel(kernel);
    }
    if(queue) {
      clReleaseCommandQueue(queue);
    }
  }

  if(program) {
    clReleaseProgram(program);
  }
  if(context) {
    clReleaseContext(context);
  }
}
#endif
