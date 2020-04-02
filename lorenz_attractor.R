library(ggplot2)
# Solve differential equations
library(deSolve)
#smoothly animate the transformation
#devtools::install_github("thomasp85/transformr")
library(transformr)
# animates ggplot2
library(gganimate)
library(wesanderson)

# Specify the wiki equations
lorenz_equations <- function(dt, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- sigma * (y - x)
    dy <- x * (rho - z) - y
    dz <- x * y - beta * z
    list(c(dx, dy, dz))
  })
}

# Solve with initial condition
lorenz_solve <- function(y, dt, params) {
  as.data.frame(
    deSolve::ode(y = y, times = dt, 
                 func = lorenz_equations, 
                 parms = params, 
                 method = "ode45")
  )
}

# Specify constants (came directly from Wiki)
constants <- c(sigma = 10, beta = 8/3, rho = 28)

# Time over which the equations are numerically solved
times <- seq(0, 25, 0.01)

# Initial conditions, in the order x, y, z
state <- c(x = 2, y = 3, z = 4)
state_next <- c(x = 2, y = 3, z = 4.1)

# Solve the diffrential equations
solution_a <- lorenz_solve(state, times, constants)
solution_b <- lorenz_solve(state_next, times, constants)
solution_a$initialization <- "(2,3,4)"
solution_b$initialization <- "(2,3,4.1)"

# Combine the two solutions into the same data frame for plotting
all_solutions <- rbind(solution_a[1:2500, ],solution_b[1:2500, ])

plot <- ggplot2::ggplot(data= all_solutions, aes(x = x, y = y, z = z, color = initialization)) +
                  theme_bw() +
                  #gg3D::axes_3D(theta = -135, phi = 14) +
                  gg3D::stat_3D(data = all_solutions[,-1], mapping = aes(x = x, y = y, z = z, color = initialization), 
                                geom = "path", size = 0.5, theta = -135, phi = 14) +
                  gg3D::stat_3D(geom = "point", theta = -135, phi = 14) +
                  coord_equal() +
                  transition_time(time) +
                  scale_colour_manual(values= wesanderson::wes_palette("Royal1", 4)) +
                  labs(title = "Lorenz Attractor Animated") +
                  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size= 25, color = "slategray4"),
                        axis.line = element_blank(), 
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(), 
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(), 
                        legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) 

# Use lower fps to slowdown particle speed on the attractor
gganimate::animate(plot, nframes = 2500, fps = 20, height= 800, width= 800) 

# Save it
gganimate::anim_save("bfe_3d_animation.gif")
