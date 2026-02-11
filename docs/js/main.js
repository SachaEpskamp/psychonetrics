// psychonetrics website — main.js

document.addEventListener('DOMContentLoaded', function () {

  // --- Mobile navigation toggle ---
  var toggle = document.querySelector('.nav-toggle');
  var menu = document.querySelector('.nav-menu');
  if (toggle && menu) {
    toggle.addEventListener('click', function () {
      menu.classList.toggle('open');
    });
    document.addEventListener('click', function (e) {
      if (!e.target.closest('.nav-container') && menu.classList.contains('open')) {
        menu.classList.remove('open');
      }
    });
  }

  // --- Highlight active nav link ---
  var page = window.location.pathname.split('/').pop() || 'index.html';
  document.querySelectorAll('.nav-link').forEach(function (link) {
    if (link.getAttribute('href') === page) {
      link.classList.add('active');
    }
  });

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
