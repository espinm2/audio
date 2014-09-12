HOMEWORK 2: CLOTH & FLUID SIMULATION

NAME:  < Max Espinoza >


TOTAL TIME SPENT:  < 20+ >
Please estimate the number of hours you spent on this assignment.


COLLABORATORS AND OTHER RESOURCES: 
List the names of everyone you talked to about this assignment
(classmates, TAs, students/instructor via LMS, etc.), and all of the
resources (books, online reference material, etc.) you consulted in
completing this assignment.

< Everyone in the Graphics Labs >

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.


DISCUSS STABILITY OF YOUR CLOTH SIMULATION:

My cloth simulation is pretty stable, it wasn't great until the implementation of the adaptive timestep.
The way I handle that is if a given particle/mass tries to move more then half the distance of a spring at rest
I reduce the time step by halve. I've had pretty good results with these. The table_colth simulation worked great,
and you even get a really cool slow down effect if you start off with a higher tempstep, matrix like. Anyways, so long
as the time step isn't large enough to rip the cloth apart my simulation is convicing. However given some test cases such as
making either the mass to light.


DESCRIBE YOUR NEW CLOTH TEST SCENE & 
THE COMMAND LINE TO RUN YOUR EXAMPLE:

I couldn't think of many large objects made of cloth, aside from a table and tent.
So i went with recreating a cloth tent, like the one at fairs as such: http://i.ytimg.com/vi/ZwfBuKktcEg/0.jpg
I think it caputred the essence of a cloth tent pretty well!

./simulation -cloth cloth_tent.txt



DISCUSS THE ACCURACY & STABILITY OF YOUR FLUID SIMULATION:

I didn't impelement the 3d velocity interpolation, or surface cells. But aside from that I think the behavoir is as expected went it comes to natural looking particles moving through fluid. I checked with the professor and we finally worked out the
last few bugs in xy interpolation. I also used the same principle I had in cloth, in which I had an adaptive timestep. Such that if given that a particle travels in one timestep a width a cell. In such a case I halve the time step, and it keeps the system from breaking. Aside from that I think it's pretty stable with the examples we were given to use.

I 


KNOWN BUGS IN YOUR CODE:
Please be concise!
Went OCD to get rid of the last few obvious ones in the fluid. The main bug I have in my code for cloth is the lack of being
able to render a diffrence between deniem and slik. They both looked the same to me when I ran them.



NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:

Include instructions for use and test cases and sample output as appropriate.
I implemented the fluid compression in more then one manner. I used the approach by the paper that used an approximation function to zero out the pressures, but I also use the method described in class. To Foster method in the paper,
on the input of a fluid file, just write in Foster before "grid" and the program will run the alternative method. I've attached a file called fluid_incompressible_foster.txt that when used should run. I like the Foster one because on average I found it was much more stable then the one used in class.

Copy & Paste Line
./simulation -fluid ../src/fluid_incompressible_foster.txt 


