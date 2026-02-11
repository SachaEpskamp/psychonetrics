// psychonetrics website — main.js

document.addEventListener('DOMContentLoaded', function () {

  // --- Mobile navigation toggle ---
  var toggle = document.querySelector('.nav-toggle');
  var menu = document.querySelector('.nav-menu');
  var dropdown = document.querySelector('.nav-dropdown');
  var dropdownToggle = document.querySelector('.dropdown-toggle');

  if (toggle && menu) {
    toggle.addEventListener('click', function () {
      menu.classList.toggle('open');
    });
  }

  // --- Click outside: close mobile menu + dropdown ---
  document.addEventListener('click', function (e) {
    if (menu && !e.target.closest('.nav-container') && menu.classList.contains('open')) {
      menu.classList.remove('open');
    }
    if (dropdown && !e.target.closest('.nav-dropdown')) {
      dropdown.classList.remove('open');
      if (dropdownToggle) dropdownToggle.setAttribute('aria-expanded', 'false');
    }
  });

  // --- Dropdown submenu ---
  if (dropdown && dropdownToggle) {
    dropdownToggle.addEventListener('click', function (e) {
      e.stopPropagation();
      var isOpen = dropdown.classList.toggle('open');
      dropdownToggle.setAttribute('aria-expanded', isOpen ? 'true' : 'false');
    });
    // Desktop: open on hover
    dropdown.addEventListener('mouseenter', function () {
      if (window.innerWidth > 768) {
        dropdown.classList.add('open');
        dropdownToggle.setAttribute('aria-expanded', 'true');
      }
    });
    dropdown.addEventListener('mouseleave', function () {
      if (window.innerWidth > 768) {
        dropdown.classList.remove('open');
        dropdownToggle.setAttribute('aria-expanded', 'false');
      }
    });
  }

  // --- Escape key closes dropdown ---
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape' && dropdown && dropdown.classList.contains('open')) {
      dropdown.classList.remove('open');
      if (dropdownToggle) {
        dropdownToggle.setAttribute('aria-expanded', 'false');
        dropdownToggle.focus();
      }
    }
  });

  // --- Highlight active nav link ---
  var page = window.location.pathname.split('/').pop() || 'index.html';
  var modelPages = ['varcov.html', 'lvm.html', 'ising.html', 'var1.html'];
  document.querySelectorAll('.nav-link').forEach(function (link) {
    if (link.getAttribute('href') === page) {
      link.classList.add('active');
    }
  });
  if (modelPages.indexOf(page) !== -1 && dropdownToggle) {
    dropdownToggle.classList.add('parent-active');
  }

  // --- Copy-code buttons ---
  document.querySelectorAll('pre[class*="language-"]').forEach(function (block) {
    var btn = document.createElement('button');
    btn.className = 'copy-btn';
    btn.textContent = 'Copy';
    btn.addEventListener('click', function () {
      var code = block.querySelector('code');
      var text = code ? code.textContent : block.textContent;
      navigator.clipboard.writeText(text).then(function () {
        btn.textContent = 'Copied!';
        btn.classList.add('copied');
        setTimeout(function () {
          btn.textContent = 'Copy';
          btn.classList.remove('copied');
        }, 2000);
      });
    });
    block.appendChild(btn);
  });

  // --- Example card expand/collapse ---
  document.querySelectorAll('.example-header').forEach(function (header) {
    header.addEventListener('click', function () {
      var body = header.nextElementSibling;
      var chevron = header.querySelector('.example-chevron');
      if (body && body.classList.contains('example-body')) {
        body.classList.toggle('open');
        if (chevron) chevron.classList.toggle('open');
      }
    });
  });

  // --- Smooth scroll for anchor links ---
  document.querySelectorAll('a[href^="#"]').forEach(function (anchor) {
    anchor.addEventListener('click', function (e) {
      var target = document.querySelector(this.getAttribute('href'));
      if (target) {
        e.preventDefault();
        target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }
    });
  });

});
