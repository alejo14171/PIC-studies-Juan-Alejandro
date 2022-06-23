#include <sys/stat.h>
// ===========================================================
//
//  Create output folders if necessary
//
// ===========================================================
void createFolders()
{
      char resultfolder[100];
      sprintf(resultfolder,"./output");
      if(mkdir (resultfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",resultfolder );
      else
        printf("--Created %s \n \n", resultfolder);

    char restfolder[100];
      sprintf(restfolder,"%s/results",resultfolder);
      if(mkdir (restfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",restfolder );
      else
        printf("--Created %s \n \n", restfolder);

     char restelfolder[100];
      sprintf(restelfolder,"%s/ion",restfolder);
      if(mkdir (restelfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",restelfolder );
      else
        printf("--Created %s \n \n", restelfolder);

      char restionfolder[100];
      sprintf(restionfolder,"%s/electron",restfolder);
      if(mkdir (restionfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",restionfolder );
      else
        printf("--Created %s \n \n", restionfolder);


      char ivdffolder[100];
      sprintf(ivdffolder,"%s/ion",resultfolder);
      if(mkdir (ivdffolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",ivdffolder );
      else
        printf("--Created %s \n \n", ivdffolder);

      char evdffolder[100];
      sprintf(evdffolder,"%s/electron",resultfolder);
      if(mkdir (evdffolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",evdffolder );
      else
        printf("--Created %s \n \n", evdffolder);
/*
      char tracefolder[100];
      sprintf(tracefolder,"%s/trace",resultfolder);
      if(mkdir (tracefolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",tracefolder );
      else
        printf("--Created %s \n \n", tracefolder);
 */
}
