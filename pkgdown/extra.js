/* Force default dark mode if no user preference is stored. */
window.addEventListener('DOMContentLoaded', function() {
  if (localStorage.getItem('theme') === null) {
    document.documentElement.setAttribute('data-bs-theme', 'dark');
    var darkButton = document.querySelector('[data-bs-theme-value="dark"]');
    if (darkButton) {
      // Simulate click to update button state and lightswitch.js internal state
      darkButton.click();
      // Remove it from localStorage so it doesn't force this preference later over system preference if the user wants auto
      localStorage.removeItem('theme');
    }
  }
});
