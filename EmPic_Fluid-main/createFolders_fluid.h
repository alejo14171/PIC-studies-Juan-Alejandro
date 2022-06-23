#include <sys/stat.h>

// ===========================================================
// 
//  Create output folders if necessary 
// 
// ===========================================================
void createFolders()
{
      char resultfolder[100];
      sprintf(resultfolder,"./results");
      if(mkdir (resultfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",resultfolder );
      else 
        printf("--Created %s \n \n", resultfolder);

      char checkfolder[100];
      sprintf(checkfolder,"./check");
      if(mkdir (checkfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",checkfolder );
      else 
        printf("--Created %s \n \n", checkfolder);


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


      char efieldsfolder[100];
      sprintf(efieldsfolder,"%s/fields",resultfolder);
      if(mkdir (efieldsfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",efieldsfolder );
      else
        printf("--Created %s \n \n", efieldsfolder);


      char seevdffolder[100];
      sprintf(seevdffolder,"%s/see",resultfolder);
      if(mkdir (seevdffolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
        printf("--Already there %s mkdir() error \n \n",seevdffolder );
      else 
        printf("--Created %s \n \n", seevdffolder);


}
