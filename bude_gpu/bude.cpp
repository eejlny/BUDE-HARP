#include <CL/cl.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <Winsock2.h>
//#include <sys/time.h>

#define POSES_PER_WI  4
#define WGSIZE        256

// OpenCL runtime configuration
cl_platform_id platform = NULL;
//unsigned num_devices = 0;
cl_device_id device; // num_devices elements
cl_context context = NULL;
cl_command_queue queue; // num_devices elements
cl_program program = NULL;
cl_kernel kernel; // num_devices elements

void listCL();
void initCL();
void releaseCL();

// High-resolution timer.
double getCurrentTimestamp() {
#ifdef _WIN32 // Windows
  // Use the high-resolution performance counter.

  static LARGE_INTEGER ticks_per_second = {};
  if(ticks_per_second.QuadPart == 0) {
    // First call - get the frequency.
    QueryPerformanceFrequency(&ticks_per_second);
  }

  LARGE_INTEGER counter;
  QueryPerformanceCounter(&counter);

  double seconds = double(counter.QuadPart) / double(ticks_per_second.QuadPart);
  return seconds;
#else         // Linux
  timespec a;
  clock_gettime(CLOCK_MONOTONIC, &a);
  return (double(a.tv_nsec) * 1.0e-9) + double(a.tv_sec);
#endif
}


void checkError(cl_int err, const char *op)
{
  if (err != CL_SUCCESS)
  {
    printf("Error during operation '%s' (%d)\n", op, err);
    releaseCL();
  }
}

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

int main(int argc, char *argv[])
{
   printf("Starting BUDE NVIDIA\n");
  cl_int err;

  //memset(&cl, 0, sizeof(cl));
  //see what is available
  listCL();

  initCL();

  cl_int natlig = sizeof(h_ligand)/sizeof(Atom);
  cl_int natpro = sizeof(h_protein)/sizeof(Atom);
  cl_int ntypes = sizeof(h_forcefield)/sizeof(FFParams);
  //cl_int nposes = sizeof(h_poses0)/sizeof(float);
  cl_int nposes = sizeof(h_poses0)/sizeof(float);
  cl_mem protein, ligand, energies, forcefield, poses[6];


  // Create buffers
  protein = clCreateBuffer(context, CL_MEM_READ_ONLY,
                           natpro*sizeof(Atom), NULL, &err);
  checkError(err, "creating protein");

  ligand = clCreateBuffer(context, CL_MEM_READ_ONLY,
                          natlig*sizeof(Atom), NULL, &err);
  checkError(err, "creating ligand");

  energies = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                            nposes*sizeof(cl_float), NULL, &err);
  checkError(err, "creating energies");

  forcefield = clCreateBuffer(context, CL_MEM_READ_ONLY,
                              ntypes*sizeof(FFParams), NULL, &err);
  checkError(err, "creating forcefield");

  for (int i = 0; i < 6; i++)
  {
    poses[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
                              nposes*sizeof(cl_float), NULL, &err);
  }

  // Write data to device
  err = clEnqueueWriteBuffer(queue, protein, CL_TRUE, 0,
                             natpro*sizeof(Atom), h_protein,
                             0, NULL, NULL);
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

  // Set kernel arguments
  err  = clSetKernelArg(kernel, 0, sizeof(cl_int), &natlig);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_int), &natpro);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &protein);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &ligand);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), poses+0);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), poses+1);
  err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), poses+2);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), poses+3);
  err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), poses+4);
  err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), poses+5);
  err |= clSetKernelArg(kernel, 10, sizeof(cl_mem), &energies);
  err |= clSetKernelArg(kernel, 11, sizeof(cl_mem), &forcefield);
  err |= clSetKernelArg(kernel, 12, ntypes*sizeof(FFParams), NULL);
  err |= clSetKernelArg(kernel, 13, sizeof(cl_int), &ntypes);
  err |= clSetKernelArg(kernel, 14, sizeof(cl_int), &nposes);
  checkError(err, "setting arguments");

  size_t global[1] = {nposes/POSES_PER_WI};
  size_t local[1] = {WGSIZE};

  double start = getCurrentTimestamp() * 1.0e+9; //ns

  // Warm-up run (not timed)
  err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
                               global, local, 0, NULL, NULL);
  checkError(err, "queuing kernel");
  err = clFinish(queue);

   double end = getCurrentTimestamp() * 1.0e+9; //ns
  checkError(err, "running kernel");

  //struct timeval tv;
  //gettimeofday(&tv, NULL);
  //double start = tv.tv_usec + tv.tv_sec*1e6;

  
  double ms = ((end-start))*1.0e-6;
  printf("Warm up time = %.1lf ms\n", ms);

  
  start = getCurrentTimestamp() * 1.0e+9; //ns


  cl_event pro_event;
  

  // Timed runs
#define RUNS 10
  for (int i = 0; i < RUNS; i++)
  {
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
                                 global, local, 0, NULL, NULL);
  }

  err = clFinish(queue);
  checkError(err, "running kernel");

   end = getCurrentTimestamp() * 1.0e+9; //ns

   

  printf("Reading results\n");

  //gettimeofday(&tv, NULL);
  //double end = tv.tv_usec + tv.tv_sec*1e6;

  // Read results
  float *results = (float *)malloc(nposes*sizeof(float));
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
  printf("GFLOP/s:         %6.3lf\n", gflops);;

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

void listCL()
{

	int i, j;
	char* value;
	size_t valueSize;
	cl_uint platformCount;
	cl_platform_id* platforms;
	cl_uint deviceCount;
	cl_device_id* devices;
	cl_uint maxComputeUnits;

	

	// get all platforms
	clGetPlatformIDs(0, NULL, &platformCount);
	//platformCount = 4;
	platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)* platformCount);
	clGetPlatformIDs(platformCount, platforms, NULL);

	//list platforms and devices

	for (i = 0; i < platformCount; i++) {
		printf(" platform number %d\n", i);
		// get all devices
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
		devices = (cl_device_id*)malloc(sizeof(cl_device_id)* deviceCount);
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

		// for each device print critical attributes
		for (j = 0; j < deviceCount; j++) {

			// print device name
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
			value = (char*)malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
			printf("%d. Device: %s\n", j + 1, value);
			free(value);

			// print hardware device version
			clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
			value = (char*)malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
			printf(" %d.%d Hardware version: %s\n", j + 1, 1, value);
			free(value);

			// print software driver version
			clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
			value = (char*)malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
			printf(" %d.%d Software version: %s\n", j + 1, 2, value);
			free(value);

			// print c version supported by compiler for device
			clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
			value = (char*)malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
			printf(" %d.%d OpenCL C version: %s\n", j + 1, 3, value);
			free(value);

			// print parallel compute units
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
				sizeof(maxComputeUnits), &maxComputeUnits, NULL);
			printf(" %d.%d Parallel compute units: %d\n", j + 1, 4, maxComputeUnits);

		}
		free(devices);

	}

	free(platforms);

	printf("\n\n\n\n");

}

void initCL()
{
#define PLATFORM 0
#define DEVICE   1
#define MAX_NUM  8

	//p0,d0 => egpu  p0,d1 => core, p1,d0 k600

  cl_int err;

  cl_platform_id platform;
  cl_platform_id platforms[MAX_NUM];
  err = clGetPlatformIDs(MAX_NUM, platforms, NULL);
  checkError(err, "getting platforms");
  platform = platforms[PLATFORM];

  cl_device_id devices[MAX_NUM];
  err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, MAX_NUM, devices, NULL);
  checkError(err, "getting devices");
  device = devices[DEVICE];

  

  char name[128];
  clGetDeviceInfo(device, CL_DEVICE_NAME, 128, name, 0);
  printf("Using device: %s\n", name);

  context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
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
  err = clBuildProgram(program, 1, &device, options, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    if (err == CL_BUILD_PROGRAM_FAILURE)
    {
      char log[16384];
      clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                            16384, log, NULL);
      printf("%s\n", log);
    }
  }
  checkError(err, "building program");

  kernel = clCreateKernel(program, "fasten_main", &err);
  checkError(err, "creating kernel");
}

#define RELEASE(func, obj) if (obj) {func(obj); obj=NULL;};
void releaseCL()
{
  RELEASE(clReleaseKernel, kernel);
  RELEASE(clReleaseProgram, program);
  RELEASE(clReleaseCommandQueue, queue);
  RELEASE(clReleaseContext, context);
}
